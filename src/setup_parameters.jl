#parent folder data
#filepath_parent = "W:\\ALL\\Theory Group\\iesh\\data\\"
filepath_parent = "/mnt/MBPC11500/braza2/data_iesh/"

#file path for initialisation files:
#filename = "C:\\Users\\Belal\\.julia\\dev\\My_iesh.jl\\data\\surface_Au111.dat"
#filename = "/home/razamaza/dev/My_iesh.jl/data/surface_Au111.dat"
#datapath = "C:\\Users\\braza2\\Documents\\GitHub\\My_iesh.jl\\data\\"
datapath = "/home/braza2/My_iesh.jl/data/"

#logging options:
#   0 -> "normal" logging
#   1 -> logging of various parameters but no logging of velocity and position vectors of bulk atoms
global const extensive_logging = 0
#conversion factors to simulation units.

#constants in SI units
global const amunit = 1.66053886e-27    #1 amu in kg
global const ev = 1.60217646e-19 # 1 eV in Joule
global const N_a = 6.02214076e23 #avogadro constant

#conversions
global const ev_kjmol = N_a * ev *1.0e-3 #ev to kj/mol
global const kjmol_seunit = 1.0e-4  #kj/mol to Dalton*angstrom^2/femtosecon^2
global const nm_to_daang = 1.0/amunit / 1.0e30 #N/m in dalton/angstrom
#physical constants
global const kb = 1.3806504e-23*1.0e-3*N_a*kjmol_seunit#	# boltzman constant [J/Kelvin] in Dalton*angstrom^2/femtosecond/Kelvin
global const hbar = 1.05457148e-34*1.0e-3*N_a*kjmol_seunit*1.0e15 #	# hbar [J*s] in Dalton*angstrom^2/femtosecond


#Simulation parameters
const global Ne = 20     #number of electrons
const global Ms = 40    #number of orbitals
const global numtraj = 100   #number of trajectories
const global tsteps = 100000
const global thop = 1   #number of timesteps between surface hops
const global twrite = 5     #number of timesteps between data writing
const global dt = 0.1  #time stepsize in Femtoseconds
#these energy parameters are used in subsequent calculations with other variables defined
# in the simulation units Angstrom, amu, femtoseconds. Therefore we need to convert
# ev to kilojoule/mole to Dalton*angstrom^2/femtosecond^2

const global e_trans_i = 0.05 * ev_kjmol * kjmol_seunit #translational energy in eV -> kilojoule/mole
const global e_vib_i = 3.1  * ev_kjmol * kjmol_seunit #vibrational energy in eV-> kilojoule/mole
const global e_rot_i = 0.00  * ev_kjmol * kjmol_seunit #translational energy in eV-> kilojoule/mole

const global r_no = 1.15077		# Initial bond length of NO
const global phi_inc = 0.0 #incident angle in radians
const global tsurf = 300.0 #temperature of surface in Kelvin
const global sim_time = tsteps * dt

#random seed used for simulation
Random.seed!(1234);

#parameters of surface etc.
const global delta_E = 7.0 * ev_kjmol * kjmol_seunit
const global sqrt_de = sqrt(delta_E)

#end simulation when scattered projectile reaches z_end
const global z_end = 11.0




#Parameters of Neutral Diabatic Matrix Element
global const A_n = 457000.052788*kjmol_seunit #Au--O: exponential repulsion: A
global const α_n = 3.75257871821 #Au--O: exponential repulsion: alpha
global const B_n = 30788.8486039*kjmol_seunit #Au--N: exponential repulsion: B
global const β_n = 2.9728905911 #Au--N: exponential repulsion: beta
global const cutoff = 10.0 #Au--N: cutoff
global const exp_beta_n_cutoff = exp(-β_n * cutoff)
global const exp_alpha_n_cutoff = exp(-α_n * cutoff)
global const F_n = 638.5*kjmol_seunit #N--O: Morse: F
global const δ_n=2.743 #N--O: Morse: delta
global const r_0_NO = 1.1507700000#N--O: Morse: cutoff

#Parameters of IONIC Diabatic Matrix Element
global const C_i = 1.25581276843  # Image potential: C
global const D_i = 347.2225355*kjmol_seunit  #Image potential: D
global const z0_i = 1.153606314 #Image potential: z0
global const A_i = 457000.052788*kjmol_seunit  #Au--O: exponential repulsion: A
global const α_i = 3.75257871821  #Au--O: exponential repulsion: alpha
global const B_i = 23.8597594272*kjmol_seunit  #Au--N: Morse: B
global const β_i=1.91014033785  #Au--N: Morse: beta
global const rN_sur_e_i = 2.38958832878 #Au--N: Morse: rN_sur_e
global const F_i = 495.9809807*kjmol_seunit   #N--O: Morse: F
global const δ_i = 2.47093477934  #N--O: Morse: delta
global const rNO_e_i = 1.29289837288  #N--O: Morse: cutoff
global const KI = 512.06425722*kjmol_seunit  #N--O: Morse: cutoff
global const exp_alpha_i_cutoff = exp(-α_i * cutoff)
global const exp_beta_i_cutoff1 = exp(-2*β_i * (cutoff - rN_sur_e_i))
global const exp_beta_i_cutoff2 = exp(-β_i * (cutoff - rN_sur_e_i))
# DEFINE PARAMETERS FOR COUPLING FUNCTION


global const coup_a_N=-70.5259924491*kjmol_seunit # Au--N: Exponential decay: a
global const coup_b_N=0.00470023958504 # Au--N: Exponential decay: b
global const coup_β_N=1.95982478112 # Au--N: Exponential decay: beta
global const coup_a_O=-16.7488672932*kjmol_seunit # Au--O: Exponential decay: a
global const coup_b_O=0.00617151653727 # Au--O: Exponential decay: b
global const coup_β_O=1.35353579356  # Au--O: Exponential decay: beta
global const Au_O_coupling_cutoff=10.0 # Au--O: Exponential decay: cutoff
global const coup_cutoff_O = coup_a_O/(1+coup_b_O*exp(coup_β_O * Au_O_coupling_cutoff))
global const coup_cutoff_N = coup_a_N/(1+coup_b_N*exp(coup_β_N * Au_O_coupling_cutoff))
#Define parameters for Au-Au interaction potential


#function to convert arrays to static/mutable staticarrays
array_to_ma(x, N, c, L) = MMatrix{N, c, Float64, L}(x)
array_to_sa(x, N, c, L) = SMatrix{N, c, Float64, L}(x)
# in Newton/meters ~ kilogram/second^2, multiply with nm_to_daang to convert to Dalton/femtosecond^2
global const α = -4.94*nm_to_daang
global const β = 17.15*nm_to_daang
global const γ = 19.40*nm_to_daang

"""u defines a transformation matrix that rotates the 111 basis vectors
([1 0 -1]/sqrt(2)],[1 -4 1]/sqrt(6),[1 1 1]/sqrt(3)) to the cartesian coordinate system
([1 0 0],[0 1 0],[0 0 1]), i.e. rotation by 45 degree around z axis and arccos(1/sqrt(3))
around x axis. """
u = permutedims([[-1.0 0.0 1.0]/sqrt(2.0);[1.0 -2.0 1.0]/sqrt(6.0);[-1.0 -1.0 -1.0]/sqrt(3.0)])
global const U_sa = array_to_sa(u, 3, 3, 9)

@doc """
    d17, d28, d39, d410, d511, d612 = def_d_matrices(α, β, γ)
dynamical matrices (force response to displacements of other atoms away from equilibrium) in cartesian basis.
The indices 17, 28, 612 refer to the 12 direct neighbors in a fcc unit cell. We define
the indices such that 1 and 7, 2 and 8 etc. have a identical dynamical matrix. For
detailed information we refer to "Thermal scattering of X-rays by crystals II.
The thermal scattering of the face-centred cubic and the close-packed hexagonal lattices
G. H. Begbie, https://doi.org/10.1098/rspa.1947.0004".
"""
function def_d_matrices(α, β, γ)
    d17 = [[α 0.0 0.0]; [0.0 β γ]; [0.0 γ β]]
    d28 = [[α 0.0 0.0]; [0.0 β -γ]; [0.0 -γ β]]
    d39 = [[β 0.0 γ]; [0.0 α 0.0]; [γ 0.0 β]]
    d410 = [[β 0.0 -γ]; [0.0 α 0.0]; [-γ 0.0 β]]
    d511 = [[β γ 0.0]; [γ β 0.0]; [0.0 0.0 α]]
    d612 = [[β -γ 0.0]; [-γ β 0.0]; [0.0 0.0 α]]
    d17 = array_to_sa(d17, 3, 3, 9)
    d28 = array_to_sa(d28, 3, 3, 9)
    d39 = array_to_sa(d39, 3, 3, 9)
    d410 = array_to_sa(d410, 3, 3, 9)
    d511 = array_to_sa(d511, 3, 3, 9)
    d612 = array_to_sa(d612, 3, 3, 9)
    return d17, d28, d39, d410, d511, d612
end
global const d17, d28, d39, d410, d511, d612 = def_d_matrices(α, β, γ)

@doc """
d1_new, d2_new, d3_new, d4_new, d5_new, d6_new = compute_d_new_basis(d17, d28, d39, d410, d511, d612)

Change of basis to 111 basis.
"""
function compute_d_new_basis(d17, d28, d39, d410, d511, d612)
    d1_new = U_sa' * d17 * U_sa
    d2_new = U_sa' * (d28* U_sa)
    d3_new = U_sa' * (d39* U_sa)
    d4_new = U_sa'* (d410* U_sa)
    d5_new = U_sa' * (d511* U_sa)
    d6_new = U_sa' * (d612* U_sa)
    return d1_new, d2_new, d3_new, d4_new, d5_new, d6_new
end
const global d1_new, d2_new, d3_new, d4_new, d5_new, d6_new = compute_d_new_basis(d17, d28, d39, d410, d511, d612)

#set masses
const global mass_arr =  SVector{3, Float64}(14.00307440, 15.99491502, 196.966548)
const global m_N = mass_arr[1]
const global m_O = mass_arr[2]
const global m_au = mass_arr[3]
const global μ = m_N*m_O/(m_N + m_O)


#create mass vector m_spread and inverse mass vector inv_m_spread
m_spread_1 =repeat(mass_arr, inner=[1, 3])
m_spread_2 = repeat([m_au], inner=[1, 3], outer=[527, 1])
const global m_spread = vcat(m_spread_1, m_spread_2)
const global inv_m_spread = 1.0 ./ m_spread

#get equilibrium position in unit cell for gold, see "Begbie 1947 https://doi.org/10.1098/rspa.1947.0004"
function get_r0()
    r0_old_basis = 0.5*a*collect([0,1,1,0,1,-1,1,0,1,-1,0,1 ,1,1,0 ,1,-1,0,0,-1,
    -1,0,-1,1,-1,0,-1,1, 0,-1,-1,-1,0,-1,1,0])
    r0_old_basis = reshape(r0_old_basis, (3,12))
    r0_new_basis = zeros(Float64, 3, 12)
    mul!(r0_new_basis, transpose(U_sa) ,r0_old_basis)
    return r0_new_basis
end

global const r0 = array_to_sa(get_r0(), 3, 12, 3*12)

@doc """
h0_temp, e_diabat_temp, vm_temp = burkey_cantrell()

Computes Continuum discretization via Burkey-Cantrell method, see
Burkey, Ronald S. / Cantrell, C. D.
Discretization in the quasi-continuum
1984-04
"""
function burkey_cantrell()
    e_diabat = zeros(Float64, Ms) #Diabat energies
    h0 = zeros(Float64, Ms+1, Ms+1) #component of diabatic hamiltonian matrix
    vm = zeros(Float64, Ms+1, Ms+1)#component of diabatic hamiltonian matrix
    gauss = zeros(Float64, Int(Ms/2), Int(Ms/2)) #recursion matrix to obtain continuum states
    for j in 1:(Int(Ms/2)-1)
        gauss[j, j+1] = Float64(j)/sqrt((2.0*Float64(j) + 1.0)*(2.0*Float64(j) - 1.0))   #equation not trivially reproducable
        gauss[j+1, j] = gauss[j, j+1]
    end

    eigvals_gauss = eigvals(gauss)
    eigvecs_gauss = eigvecs(gauss)

    for j in 2:(Int(Ms/2) + 1)
        h0[j, j] = delta_E/4.0*eigvals_gauss[j-1] - delta_E/4.0     #???
         h0[j+Int(Ms/2), j+Int(Ms/2)] = delta_E/4.0*eigvals_gauss[j-1] + delta_E/4.0
        e_diabat[j-1] = h0[j, j]
        e_diabat[j + Int(Ms/2) - 1] = h0[j+Int(Ms/2), j+Int(Ms/2)]
    end

    vm[1, 2:Int(Ms/2) + 1] = @. sqrt(delta_E/4.0 * (2.0*eigvecs_gauss[1, :]^2))  #??? why sqrt? where is weight function?
    vm[2:Int(Ms/2) + 1, 1] = @. sqrt(delta_E/4.0 * (2.0*eigvecs_gauss[1, :]^2))
    vm[1, (2 + Int(Ms/2)):(Ms + 1)] = @. sqrt(delta_E/4.0 * (2.0*eigvecs_gauss[1, :]^2))
    vm[(2 + Int(Ms/2)):(Ms + 1), 1] = @. sqrt(delta_E/4.0 * (2.0*eigvecs_gauss[1, :]^2))

    return h0, e_diabat, vm
end

h0_temp, e_diabat_temp, vm_temp = burkey_cantrell()
global const h0 = h0_temp
global const e_diabat = SVector{Ms, Float64}(e_diabat_temp)
global const vm = vm_temp

#auxilliary logical functions
leaq(a,b) = (a <= b) || (a ≈ b)
geaq(a, b) = (a >= b) || (a ≈ b)
