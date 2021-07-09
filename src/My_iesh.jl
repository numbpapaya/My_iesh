module My_iesh


curr_vers = "v.0.85"
using LinearAlgebra
using Random
using StaticArrays
using Plots
using HDF5
using Statistics
using Dates

include("read_au_data.jl");
include("setup_parameters.jl");
include("types.jl");
include("misc_func.jl");
include("constructor.jl");
include("get_neighbors.jl");
include("get_energies.jl");
include("get_forces.jl");
include("propagation.jl");
include("run_trajectory.jl");


s = simulation_init()
simulation_constructor_x_v!(s)
simulation_constructor_nn!(s)
simulation_constructor_x_300K!(s)
simulation_constructor_energy(s)
simulation_constructor_force(s)
# propagate_init!(s)
# simulate!(s)
s.H .= Symmetric(h0 .+ vm .* s.hp[3] ./ sqrt_de) #check
s.H[1, 1] = s.H[1, 1] + s.hp[2] - s.hp[1] #check
s.λ .= eigvals(s.H)
s.Γ .= eigvecs(s.H)#check seems okay
# s.Γ2 .= deepcopy(s.Γ)
# s.Γ2[:, 1] .= -s.Γ2[:, 1]
# s.Γ2[:, 2] .= -s.Γ2[:, 2]
# s.Γ2[:, 4] .= -s.Γ2[:, 4]
# s.Γ2[:, 1] .= -s.Γ2[:, 1]
# s.Γ2[:, 10] .= -s.Γ2[:, 10]
for j in 1:Ms+1
    s.dhdea[:, j] .= view(s.Γ, 1, :) .* s.Γ[1, j]
end
temp = vm * s.Γ #needs preallocation
mul!(s.dhdv, transpose(s.Γ), temp)
s.dhdv .= s.dhdv ./ sqrt_de


generate_therm_surf!(s)

#set initial wavefunction
s.ψ .= convert(Array{ComplexF64, 2}, view(s.Γ, :, s.surfp)) #BUUUUUUUG

# KE = sum(0.5 ./m_spread .* s.v.^2)
# PE = s.hp[1] + sum(s.Γ[s.surfp])

#calculate initial forces
for j in 1:Ne
    s.F .= s.F .+ (s.dhdea[s.surfp[j], s.surfp[j]] .* (s.dhp_ion .- s.dhp_neutral) .+
    s.dhdv[s.surfp[j], s.surfp[j]] .* s.dhp_coup)

    # @views s.F .= s.F .+ s.dhdv[s.surfp[j], s.surfp[j]] .* s.dhp_coup
end
s.F .= s.F .+ s.dhp_neutral
#fix position of the last layer of atoms
s.F[399:530, :] .= 0.0
s.v[399:530, :] .= 0.0

#convert to Forces
s.F .= -s.F

n = 2

@inbounds for j in 1:Ne
    s.ϕ[:, j] = view(s.Γ, :, s.surfp[j])
end


#Calculate nonadiabatic coupling matrix DM between adiabatic orbitals
# Calculate nonadiabatic coupling elements DM between adiabatic orbitals s.dm
# phipsi = <phi_k|psi> is overlap between "current" adiabatic state |phi_k> and electronic state |psi>
temp_1 = transpose(s.ϕ) * s.ψ   #preallocation
s.phipsi[1] = det(temp_1)
# Calculate akk as defined in Tully JCP 1990
s.akk[1] = abs(s.phipsi[1])^2 #corrected error

# Calculate nonadiabatic coupling elements DM between adiabatic orbitals

temp_2 = s.dhp_ion .- s.dhp_neutral #preallocation
temp_2 .= temp_2 .* s.v
temp_2_sum = sum(temp_2)
s.dm .= temp_2_sum .* s.dhdea

temp_3 = s.dhp_coup .* s.v
temp_3_sum = sum(temp_3)
s.dm .= s.dm .+ temp_3_sum .* s.dhdv

temp_rep1 = repeat(s.λ, 1, Ms+1)
temp_rep2 = temp_rep1' .- temp_rep1
s.dm .= s.dm ./ temp_rep2
for j in 1:Ms+1
    s.dm[j,j] = zero(eltype(s.dm))
end #should be ok
# Get occupied (particle), unoccupied (hole) states corresponding to current surface specified by surfp
# surfp is Ne x 1 array of particle states
# surfh is (Ms+1-Ne) x 1 array of hole states
s.occnum .= 0
for j in 1:Ne
    s.occnum[s.surfp[j]] = 1
end
k = 1
for j in 1:Ms+1
    if s.occnum[j] == 0
        s.surfh[k] = j
        k = k+1
    end
end #should be ok
# Get b_{kl} where k is current many-electron adiabatic state and l is possible new
# 	many-electron adiabatic state
# k corresponds to occupation of states specified by surfp
# l corresponds to occupation of states specified by a one electron-hole-pair excitation out of surfp
rtemp = abs(s.phipsi[1].^2)
if geaq(rtemp, 1.0)
    rtemp = 1.0
end
for jp in 1:Ne
    for jh in 1:(Ms+1-Ne)
        s.blk2[jp, jh] = 2.0 * abs(s.phipsi[1] * s.dm[s.surfp[jp], s.surfh[jh]])
    end
end
#
s.Pbmaxest[1] = sqrt(sum(s.blk2.^2)) * sqrt(1.0 - rtemp)/s.akk[1] * dt * Float64(thop) #should be ok
# Calculate upper bound of maximum hopping probability, Pbmaxest
#beginn hopping
#Generate random number for surface hopping
hoprand = 1e-25
#! Only attempt hop if hoprand < Pbmaxest to avoid expensive calculation
#! of real hopping elements blk
# if hoprand < s.Pbmaxest[1]
#     hopping!(s, hoprand, n)
# end
for jp in 1:Ne
    for jh in 1:Ms+1-Ne
        s.surfpnew .= s.surfp
        s.surfpnew[jp] = s.surfh[jh]
        for j in 1:Ne
            s.ϕ[:, j] = view(s.Γ, :, s.surfpnew[j])
        end
        #Ctemp = <phi_l|psi> is overlap between "new" adiabatic state
        #|phi_l> and electronic state |psi>
        mul!(s.ctemp1, transpose(s.ϕ) , s.ψ)
        ctemp = det(s.ctemp1)
        s.akl[1] = s.phipsi[1] * conj(ctemp)

        s.blk[jp, jh] = 2.0 * real(s.akl[1] * s.dm[s.surfp[jp], s.surfh[jh]])
    end
end

rtemp = 0.0
for jh in 1:Ms+1-Ne
    for jp in 1:Ne
        s.Pb[(jh-1) * Ne + jp] = rtemp
        s.Pb[(jh-1) * Ne + jp] = s.Pb[(jh-1)*Ne + jp] + s.blk[jp, jh]*(s.blk[jp, jh] > 0.0)
        rtemp = s.Pb[(jh - 1) * Ne + jp]
    end
end
s.Pb .= s.Pb ./ s.akk[1] * dt * Float64(thop)
# s.Pb .= collect([0.000000000000000,  0.000000000000000,  0.000000000000000,
#  8.478212327683901E-021,  1.687455563337645E-020,  1.687455563337645E-020,
#  1.687455563337645E-020,  1.687455563337645E-020,  2.589520627005461E-020,
#  2.589520627005461E-020,  3.289326919342432E-020,  3.289326919342432E-020,
#  3.572116188139076E-020,  3.637427690499731E-020,  5.227210776662324E-020,
#  5.993517139771526E-020,  6.111176661724437E-020,  6.111176661724437E-020,
#  7.700507786302974E-020,  9.041006730408093E-020,  9.041006730408093E-020,
#  9.127461762941845E-020,  9.127461762941845E-020,  9.127461762941845E-020,
#  9.127461762941845E-020,  9.128878344050694E-020,  9.136101176603941E-020,
#  9.142436034399336E-020,  9.242435982197457E-020,  9.253046132491244E-020])
# s.Pb .= collect([0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000,
# 2.095224393686503E-020,  2.095224393686503E-020,  2.095224393686503E-020,
# 2.095224393686503E-020,  4.354271952221061E-020,  9.825025048380536E-020,
# 9.825025048380536E-020,  9.825025048380536E-020,  9.825025048380536E-020,
# 9.825025048380536E-020,  9.825025048380536E-020,  1.105728131822829E-019,
# 1.105728131822829E-019,  1.275683194038428E-019,  1.384812974961126E-019,
# 1.580588791119147E-019,  1.634859067918526E-019,  1.890152151354714E-019,
# 1.938539975964183E-019,  2.033166575331342E-019,  2.033166575331342E-019,
# 2.033166575331342E-019,  2.033166575331342E-019,  2.054171378029279E-019,
# 2.054171378029279E-019,  2.054171378029279E-019,  2.162098475751005E-019])
# s.Pb .= collect([0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000,
#  3.099216181412439E-019,  3.099216181412439E-019,  3.099216181412439E-019,
#  3.099216181412439E-019,  6.432578889648027E-019,  1.454056158191952E-018,
#  1.454056158191952E-018,  1.454056158191952E-018,  1.454056158191952E-018,
#  1.454056158191952E-018, 1.454056158191952E-018,  1.638061187117117E-018,
#  1.638061187117117E-018,  1.889906289426600E-018,  2.052251515744369E-018,
#  2.344834057520700E-018,  2.426252809278021E-018,  2.804199460769983E-018,
#  2.876176861605881E-018,  3.017497247784225E-018,  3.017497247784225E-018,
#  3.017497247784225E-018,  3.017497247784225E-018,  3.048933794061390E-018,
#  3.048933794061390E-018,  3.048933794061390E-018,  3.212577971252680E-018])
s.Pb .= collect([0.000000000000000E+000,  0.000000000000000E+000,  0.000000000000000E+000,
6.195316619316498E-019,  6.195316619316498E-019,  6.195316619316498E-019,
6.195316619316498E-019,  1.285809856224723E-018,  2.906705345022110E-018,
2.906705345022110E-018,  2.906705345022110E-018,  2.906705345022110E-018,
2.906705345022110E-018,  2.906705345022110E-018,  3.274655873475869E-018,
3.274655873475869E-018,  3.778125543022411E-018,  4.102719606724644E-018,
4.687807462958685E-018,  4.850645919084630E-018,  5.606182581092673E-018,
5.750094215901735E-018,  6.032690737338966E-018,  6.032690737338966E-018,
6.032690737338966E-018,  6.032690737338966E-018,  6.095558782909743E-018,
6.095558782909743E-018,  6.095558782909743E-018,  6.422972389478660E-018])


if hoprand < s.Pb[Ne*(Ms+1-Ne)]
    hopstate = one(Int64)
    while s.Pb[hopstate] <= hoprand
        hopstate = hopstate + one(Int64)
    end
    jh = Int(floor((hopstate - 1) / Ne)) + one(Int64)
    jp = ((hopstate - one(Int64)) % Ne) + one(Int64)

    @views temp = s.dhdea[s.surfp[jp], s.surfh[jh]] .* (s.dhp_ion .- s.dhp_neutral)
    @views temp .= temp .+ s.dhdv[s.surfp[jp], s.surfh[jh]] .* s.dhp_coup
    temp2 = sqrt(sum(temp.^2)/sum(s.v .^2))
    temp .= temp ./ temp2
    temp3 = sum(temp .* s.v)
    temp .= temp .* copysign(1.0, temp3)
    rhop = deepcopy(temp)


    @views temp = s.dhdea[s.surfp[jp], s.surfh[jh]] .* (s.dhp_ion .- s.dhp_neutral)
    @views temp .= temp .+ s.dhdv[s.surfp[jp], s.surfh[jh]] .* s.dhp_coup
    temp2 = sqrt(sum(temp.^2)/sum(s.v .^2))
    temp .= temp ./ temp2
    temp3 = sum(temp .* s.v)
    temp .= temp .* copysign(1.0, temp3)
    rhop = deepcopy(temp)




    s.surfpnew .= s.surfp
    s.surfpnew[jp] = s.surfh[jh]
    e_el_old = sum(s.λ[s.surfp])
    e_el_new = sum(s.λ[s.surfpnew])
    k_old = 0.5*sum(m_spread .* s.v .^2)
    k_rhop = 0.5 * sum(m_spread .* rhop .^2)
    k_temp = sum(m_spread .* s.v .* rhop)
    bool_e = k_temp^2 - 4.0*k_rhop * (e_el_new - e_el_old) > 0.0


    s.attnum .= s.attnum .+ one(Int64)
    if bool_e == true
        rtemp = s.vscale[1]
        s.vscale[1] = (k_temp - sqrt(k_temp^2 - 4.0 * k_rhop * (e_el_new - e_el_old)))/2.0/k_rhop
        @views k_no = 0.50 * sum(view(m_spread, 1:2, :).* view(s.v, 1:2, :).^2)
        @views k_au = 0.50 * sum(view(m_spread, 3:N+2, :) .* view(s.v, 3:N+2, :).^2)
        s.v .= s.v .- s.vscale[1] .* rhop
        s.storage_deltaKNO[s.exnum[1]+1] = k_no - 0.5*sum(view(m_spread, 1:2, :) .* view(s.v, 1:2, :) .^2)
        s.storage_deltaKAu[s.exnum[1]+1] = k_au - 0.5*sum(view(m_spread, 3:N+2, :) .* view(s.v, 3:N+2, :) .^2)
        s.storage_state[2*s.exnum[1] + 1] = s.surfp[jp]
        s.storage_state[2*s.exnum[1] + 2] = s.surfh[jh]
        s.surfp .= s.surfpnew
        s.storage_hoptimes[s.exnum[1]+1] = n
        s.exnum .= s.exnum .+ one(Int64)
    end
end
uu = MVector{Ms+1, ComplexF64}(zeros(ComplexF64, Ms+1))
uuu = zeros(ComplexF64, Ms+1, Ne)
s.ψ .= transpose(s.Γ) * s.ψ

uu .= exp.(-1im * s.λ * dt/hbar)
uuu .= repeat(uu, 1, Ne)
s.ψ .= uuu .* s.ψ
s.ψ .= s.Γ * s.ψ


s.vdot .= inv_m_spread .* s.F


s.x .= s.x .+ s.v * dt .+ 0.50*s.vdot * dt^2

s.Δ_no .=  s.x[1, :] - s.x[2, :]


s.vtemp .= s.v .+ 0.50*s.vdot*dt #ARRAY ALLOCATION !!!!

s.hp .= get_E_all(s)    #update energy surface
get_F_all(s)            #update partial derivatives dhp_i
s.H .= h0 .+ vm .* s.hp[3] ./sqrt_de
s.H[1, 1] = s.H[1, 1] + s.hp[2] - s.hp[1]
Γ_hold = deepcopy(s.Γ) #ARRAY ALLOCATION !!!! solve later
s.λ .= eigvals(s.H)
s.Γ .= eigvecs(s.H)
temp = zeros(Float64, Ms+1)
@inbounds for j in 1:Ms+1
    @views temp .= abs.(transpose(Γ_hold) * s.Γ[:, j])
    Γmaxloc = argmax(temp)
    temp2 = dot(view(Γ_hold,:, Γmaxloc), view(s.Γ, :, j))
    s.Γ[:, j] = s.Γ[:, j]/copysign(1.0, temp2)
end

s.blk .= 0.0
@inbounds for j in 1:Ms+1
    s.dhdea[:, j] = view(s.Γ, 1, :) .* s.Γ[1, j]
end
temp = (vm * s.Γ) #ARRAY ALLOCATION !!!! solve later
mul!(s.dhdv, transpose(s.Γ), temp)
s.dhdv .= s.dhdv ./sqrt_de

s.F .= 0
update_forces!(s)
s.F[399:530, :] .= 0
s.v[399:530, :] .= 0
s.F .= -s.F

s.v .= s.vtemp .+ 0.5 *inv_m_spread .* s.F * dt






end























# Propagate electronic Hamiltonian
uu = MVector{Ms+1, ComplexF64}(zeros(ComplexF64, Ms+1))
uuu = zeros(ComplexF64, Ms+1, Ne)
propagate_hamiltonian!(s, uu, uuu)
#storage
if logopt == 1
    storage_sim1!(s, n)
end
#propagate nuclear motion
propagate_nuclear!(s)
#get adiabatic eigenvalues, eigenvectors
get_eigen_sim!(s)
# Get dH / d Ep matrices where Ep is parameter in Newns-Anderson Hamiltonian
# (i.e. Ep = V or Ep = E_a = (EI-EN) )
get_dhdea_dhdv_loop!(s)

# Calculate forces
s.F .= 0
update_forces!(s)
s.F[399:530, :] .= 0
s.v[399:530, :] .= 0
s.F .= -s.F

















end
