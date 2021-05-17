#module My_iesh

using LinearAlgebra
using Distances
using PyCall
using Random
using Profile
using ProfileView
using StaticArrays
#setup and definitions
include("setup_parameters.jl")
include("misc_functions.jl")

#begin reading data
include("read_au_data.jl")
#find neighbors in lattice
include("neighbors.jl") #needs testing!!!!

#define functions for energy and force calculations, yielding diabatic PES and define_forces
include("define_energies.jl")
include("define_forces.jl")

#initialize continuum representation
include("continuum_discretization.jl")
include("time_loop.jl")
#further initializations
include("initialize.jl")

#equilibrium position of 12 next neighbors in lattice
global const r0 = get_r0()
#compute continuum discretization
h0_temp, e_diabat_temp, vm_temp = burkey_cantrell()
global const h0 = h0_temp
global const e_diabat = e_diabat_temp
global const vm = vm_temp



for traj in 1:numtraj
    #begin main loop over trajectories
    #initialize arrays
    storage_x = zeros(Float64, 6, tsteps)
    storage_v = zeros(Float64, 6, tsteps)
    storage_e = zeros(Float64, 5, tsteps) #store data electronic excitation energy
    storage_aop = zeros(Float64, tsteps) #store data adsorbate orbital populations
    storage_op = zeros(Float64, Ms+1, tsteps) #store data orbital populations
    storage_xno = zeros(Float64, 6, tsteps) #store NO position
    storage_vno = zeros(Float64, 6, tsteps) #store NO velocity
    storage_phop = zeros(Float64, tsteps) #store hopping probability
    storage_e_fric = zeros(Float64, tsteps) #Stores <E_f|v(t).d|E_f>
    storage_langevin = zeros(Float64, 6, tsteps) #stores comps of the Langevin Force
    #Surface occupation arrays
    surfp = zeros(Int, Ne)
    surfpinit = zeros(Int, Ne)
    surfh = zeros(Int, Ms+1-Ne)
    surfpnew = zeros(Int, Ne)
    occnum =  zeros(Int, Ms+1)

    dhp = zeros(Float64, 3, 3 * (N+2)) #Hamiltonian parameter derivatives
    dhdea = zeros(Float64, Ms+1, Ms+1) #dh/dEa
    F = zeros(Float64, N+2, 3) #forces on atoms
    #Hopping parameters
    blk = zeros(Float64, Ne, Ms+1-Ne)
    blk2 = zeros(Float64, Ne, Ms+1-Ne)
    Pb = zeros(Float64, Ne*(Ms+1-Ne))

    #m_spread for vectorized stuff
    m_spread_1 =repeat(mass_arr, inner=[1, 3])
    m_spread_2 = repeat([m_au], inner=[1, 3], outer=[527, 1])
    m_spread = vcat(m_spread_1, m_spread_2)

    trajzmin = 10.0
    trajtheta = 0.0

    psi = zeros(ComplexF64, Ms+1, Ne) #wavefunction
    phi = zeros(ComplexF64, Ms+1, Ne) #adiabatic states


    #initialize
    x, v = initialize()
    nn = nn_lattice(x)

    hp = get_E_all(x) #Hamiltonian parameters
    dhp = get_F_all(x)

    H = h0 + vm * hp[3] / sqrt(delta_E)
    H[1, 1] = H[1, 1] + hp[2] - hp[1]
    eigval_H = eigvals(H)
    eigvec_H = eigvecs(H)

    for i in 1:Ms+1
        dhdea[:, i] = @. eigvec_H[1, :] * eigvec_H[1, i]
    end
    dhdv = eigvec_H' * (vm * eigvec_H)/sqrt(delta_E)

    surfp, occnum, surfh, surfpinit = generate_therm_surf!(surfp, occnum, surfh, surfpinit)

    #set initial wavefunction
    psi =  psi + eigvec_H[:, surfp] #add psi to make psi complex

    KE = 0.5 * (
        1/mass_arr[1] * sum(v[1,:].^2) + 1/mass_arr[2] * sum(v[2,:].^2) +
        1/mass_arr[3] * sum(v[3:end,:].^2)
     )
    PE = hp[1] + sum(eigval_H[surfp])

    #calculate initial forces
    for j in 1:Ne
        temp_arr = @. dhdea[surfp[j], surfp[j]] * (dhp[2, :] - dhp[1, :])
        temp_arr2 = reshape(temp_arr, (3, N+2))
        temp_arr3 = permutedims(temp_arr2, (2, 1))
        F = F + temp_arr3

        temp_arr = @. dhdv[surfp[j], surfp[j]] * dhp[3, :]
        temp_arr2 = reshape(temp_arr, (3, N+2))
        temp_arr3 = permutedims(temp_arr2, (2, 1))
        F = F + temp_arr3
    end
    temp_arr = reshape(dhp[1,:], (3, N+2))
    temp_arr2 = permutedims(temp_arr, (2, 1))
    F = F + temp_arr2

    #fix last layer of grid
    F[399:end, :] .= 0.0
    F = -F

    rtemp = sqrt(sum((x[1, :] - x[2, :]).^2))
    storage_e[1, 1] = 0.5 * (mass_arr[1] + mass_arr[2]) *
     sum(((mass_arr[1].*v[1, :] + mass_arr[2].*v[2, :])/(mass_arr[1] + mass_arr[2])).^2)
    T_vib = 0.50 * μ * sum(((v[1, :] - v[2, :]).*(x[1, :] - x[2, :])/rtemp).^2)
    U_vib = F_n *(1.0 - exp(-δ_n * (rtemp -  r_0_NO)))^2
    E_vib = U_vib + T_vib
    storage_e[2, 1] = E_vib
    T_tot = 0.5 * mass_arr[1] * sum(v[1, :].^2) + 0.5*mass_arr[2]*sum(v[2, :].^2)
    storage_e[3, 1] = T_tot - T_vib - storage_e[1, 1]
    storage_e[4, 1] = sum(eigval_H[surfp])

    storage_x[1:3, 1] = x[1, :]
    storage_x[4:6, 1] = x[2, :]
    storage_v[1:3, 1] = v[1, :]
    storage_v[4:6, 1] = v[2, :]

    for n in 2:tsteps
        time_loop(n)
    end #time loop
    println("finished simulation")
end # end trajectory loop





#ProfileView.view()

#end
