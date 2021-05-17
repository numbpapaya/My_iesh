#Seed

Random.seed!(1234);

#iesh or electronic friction, set 0 for iesh and 1 for electronic friction
global const latticeopt = 0

#conversion factors to simulation units.

global const conv1 = 6.02214076 # Newton/meter * Angstrom^2 in kilojoule/mole
global const conv2 = 96.48 # eV in kilojoule per mole
global const conv3 = 1e-4 #kilojoule per mole in amu*angstrom^2/femtosecond^2
#constants in SI units
global const Å = 1.0e-10 # 1 Angström in meter
global const ev = 1.60217646e-19 # 1 eV in Joule
global const amu = 1.66053886e-27	# 1 amu in kg
global const hbar = 6.582119569e-16*conv2*conv3*1e15	# hbar in amu*angstrom^2/femtosecond^2*fs
global const fs = 1.0e-15  	# fs in s

global const kjpermol = 1.660538863e-21 # 1 kJ/mol in Joules
global const N_a = 602214076000000000000000


global const kb = 0.08617*1e-3 * conv2 * conv3 	# boltzman constant [mev/Kelvin] in simulation units []

#Simulation parameters
const global Ne = 20     #number of electrons
const global Ms = 40     #number of orbitals
const global numtraj = 1   #number of trajectories
const global tsteps = 1000
const global thop = 1   #number of timesteps between surface hops
const global twrite = 5     #number of timesteps between data writing
const global dt = 0.5  #time stepsize in Femtoseconds
#these energy parameters are used in subsequent calculations with other variables defined
# in the simulation units Angstrom, amu, femtoseconds. Therefore we need to convert
# kilojoule/mole to amu*angstrom^2/femtosecond^2
const global r_no = 1.15077		# Initial bond length of NO
const global e_trans_i = 0.05*conv2*conv3 #translational energy in eV -> kilojoule/mole
const global e_vib_i = 3.1*conv2*conv3 #vibrational energy in eV-> kilojoule/mole
const global e_rot_i = 0.05*conv2*conv3 #translational energy in eV-> kilojoule/mole
const global phi_inc = 0.0 #incident angle in radians
const global tsurf = 300.0 #temperature of surface in Kelvin
const global sim_time = tsteps * dt



#parameters of surface etc.
const global delta_E = 7.0 *conv2*conv3
const global z_end = 11.0




#Parameters of Neutral Diabatic Matrix Element
global const A_n = 457000.052788*conv3 #Au--O: exponential repulsion: A
global const α_n = 3.75257871821 #Au--O: exponential repulsion: alpha
global const B_n = 30788.8486039*conv3 #Au--N: exponential repulsion: B
global const β_n = 2.9728905911 #Au--N: exponential repulsion: beta
global const cutoff = 10.0 #Au--N: cutoff
global const exp_beta_n_cutoff = exp(-β_n * cutoff)
global const exp_alpha_n_cutoff = exp(-α_n * cutoff)
global const F_n = 638.5000000000*conv3 #N--O: Morse: F
global const δ_n=2.7430000000 #N--O: Morse: delta
global const r_0_NO = 1.1507700000#N--O: Morse: cutoff

#Parameters of IONIC Diabatic Matrix Element

global const C_i = 1.25581276843  # Image potential: C
global const D_i = 347.2225355*conv3  #Image potential: D
global const z0_i = 1.153606314 #Image potential: z0
global const A_i = 457000.052788*conv3  #Au--O: exponential repulsion: A
global const α_i = 3.75257871821  #Au--O: exponential repulsion: alpha
global const B_i = 23.8597594272*conv3  #Au--N: Morse: B
global const β_i=1.91014033785  #Au--N: Morse: beta
global const rN_sur_e_i = 2.38958832878 #Au--N: Morse: rN_sur_e
global const F_i = 495.9809807*conv3   #N--O: Morse: F
global const δ_i = 2.47093477934  #N--O: Morse: delta
global const rNO_e_i = 1.29289837288  #N--O: Morse: cutoff
global const KI = 512.06425722*conv3  #N--O: Morse: cutoff
global const exp_alpha_i_cutoff = exp(-α_i * cutoff)
global const exp_beta_i_cutoff1 = exp(-2*β_i * (cutoff - rN_sur_e_i))
global const exp_beta_i_cutoff2 = exp(-β_i * (cutoff - rN_sur_e_i))
# DEFINE PARAMETERS FOR COUPLING FUNCTION


global const coup_a_N=-70.5259924491*conv3 # Au--N: Exponential decay: a
global const coup_b_N=0.00470023958504 # Au--N: Exponential decay: b
global const coup_β_N=1.95982478112 # Au--N: Exponential decay: beta
global const coup_a_O=-16.7488672932*conv3 # Au--O: Exponential decay: a
global const coup_b_O=0.00617151653727 # Au--O: Exponential decay: b
global const coup_β_O=1.35353579356  # Au--O: Exponential decay: beta
global const Au_O_coupling_cutoff=10.0 # Au--O: Exponential decay: cutoff
global const coup_cutoff_O = coup_a_O/(1+coup_b_O*exp(coup_β_O * Au_O_coupling_cutoff))
global const coup_cutoff_N = coup_a_N/(1+coup_b_N*exp(coup_β_N * Au_O_coupling_cutoff))
#Define parameters for Au-Au interaction potential


# in Newton/meters, multiply with conv1 so that resulting units of energy will
# be in kilojoule per mole
global const α = -4.94*conv1
global const β = 17.15*conv1
global const γ = 19.40*conv1

function def_d_matrices(α, β, γ)
    #see paper from 1947 for definitions
    d17 = [[α 0.0 0.0]; [0.0 β γ]; [0.0 γ β]]
    d28 = [[α 0.0 0.0]; [0.0 β -γ]; [0.0 -γ β]]
    d39 = [[β 0.0 γ]; [0.0 α 0.0]; [γ 0.0 β]]
    d410 = [[β 0.0 -γ]; [0.0 α 0.0]; [-γ 0.0 β]]
    d511 = [[β γ 0.0]; [γ β 0.0]; [0.0 0.0 α]]
    d612 = [[β -γ 0.0]; [-γ β 0.0]; [0.0 0.0 α]]

    return d17, d28, d39, d410, d511, d612
end
const global d17, d28, d39, d410, d511, d612 = def_d_matrices(α, β, γ)

function compute_d_new_basis(d17, d28, d39, d410, d511, d612)
    U = permutedims([[-1.0 0.0 1.0]/sqrt(2.0);[1.0 -2.0 1.0]/sqrt(6.0);[-1.0 -1.0 -1.0]/sqrt(3.0)])
    d1_new = U' * d17
    d2_new = U' * d28
    d3_new = U' * d39
    d4_new = U' * d410
    d5_new = U' * d511
    d6_new = U' * d612
    return d1_new, d2_new, d3_new, d4_new, d5_new, d6_new
end
const global d1_new, d2_new, d3_new, d4_new, d5_new, d6_new = compute_d_new_basis(d17, d28, d39, d410, d511, d612)

#set masses
const global mass_arr =  [14.00307440, 15.99491502, 196.966548]
const global m_N = mass_arr[1]
const global m_O = mass_arr[2]
const global m_au = mass_arr[3]
const global μ = m_N*m_O/(m_N + m_O)
