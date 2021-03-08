#Parameters of Neutral Diabatic Matrix Element
global const A_n = 457000.052788*kjpermol #Au--O: exponential repulsion: A
global const α_n = 3.75257871821/Å #Au--O: exponential repulsion: alpha
global const Au_O_cutoff_n= 10.0 * Å #Au--O: cutoff
global const B_n = 30788.8486039d0*kjpermol #Au--N: exponential repulsion: B
global const β_n = 2.9728905911/Å #Au--N: exponential repulsion: beta
global const Au_N_cutoff_n = 10.0  * Å #Au--N: cutoff
global const F_n = 638.5000000000*kjpermol #N--O: Morse: F
global const δ_n=2.7430000000/Å #N--O: Morse: delta
global const r_0_NO = 1.1507700000*Å #N--O: Morse: cutoff

#Parameters of Neutral Diabatic Matrix Element

global const C_i = 1.25581276843 * Å  # Image potential: C
global const D_i = 347.2225355 * kjpermol * Å  #Image potential: D
global const z0_i = 1.153606314 * Å  #Image potential: z0
global const A_i = 457000.052788 * kjpermol  #Au--O: exponential repulsion: A
global const α_i = 3.75257871821 / Å  #Au--O: exponential repulsion: alpha
global const Au_O_cutoff_i = 10 / Å #Au--O: cutoff
global const B_i = 23.8597594272 * kjpermol  #Au--N: Morse: B
global const β_i=1.91014033785 / Å  #Au--N: Morse: beta
global const rN_sur_e_i = 2.38958832878 * Å  #Au--N: Morse: rN_sur_e
global const Au_N_cutoff_i = 10 * Å  #Au--N: cutoff
global const F_i = 495.9809807 * kjpermol   #N--O: Morse: F
global const δ_i = 2.47093477934 / Å  #N--O: Morse: delta
global const rNO_e_i = 1.29289837288 / Å  #N--O: Morse: cutoff
global const KI = 512.06425722 * kjpermol  #N--O: Morse: cutoff

# DEFINE PARAMETERS FOR COUPLING FUNCTION


global const coup_a_N=-70.5259924491 * kjpermol  # Au--N: Exponential decay: a
global const coup_b_N=0.00470023958504 # Au--N: Exponential decay: b
global const coup_β_N=1.95982478112 / Å # Au--N: Exponential decay: beta
global const Au_N_coupling_cutoff=10 * Å # Au--N: Exponential decay: cutoff
global const coup_a_O=-16.7488672932 * kjpermol # Au--O: Exponential decay: a
global const coup_b_O=0.00617151653727 # Au--O: Exponential decay: b
global const coup_β_O=1.35353579356 / Å  # Au--O: Exponential decay: beta
global const Au_O_coupling_cutoff=10 * Å # Au--O: Exponential decay: cutoff

#Define parameters for Au-Au interaction potential
global const α = -4.94
global const β = 17.15
global const γ = 19.40
