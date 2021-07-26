@doc """
    simulation_init()

    Constructs the simulation type with zeros or undefined values
"""
function simulation_init()
    A = Ne
    C = Ms+1
    D = Ms-Ne+1
    if extensive_logging == true
    return Simulation{A, C, D}(
    zeros(Float64, N+2, 3), #x
    zeros(Float64, N+2, 3), #v
    MVector{3, Float64}(0.0, 0.0, 0.0),    #Δ_no
    zeros(Int64, N, 12), #nn_arr
    zeros(Float64, N+2, 3), #dhp_neutral
    zeros(Float64, N+2, 3), #dhp_ion
    zeros(Float64, N+2, 3), #dhp_coup
    MVector{3, Float64}(0.0, 0.0, 0.0), #hp
    zeros(Float64, N+2, 3), #F
    MVector{A, Int16}(zeros(Int16, Ne)), #surfp
    MVector{C, Int16}(zeros(Int16, Ms+1)), #occnum
    MVector{D, Int16}(zeros(Int16, Ms-Ne+1)), #surfh
    MVector{A, Int16}(zeros(Int16, Ne)), #surfpinit
    MVector{A, Int16}(zeros(Int16, Ne)), #surfpnew
    MVector{1, Float64}(10.0), # trajzmin
    MVector{1, Float64}(0.0), #trajtheta
    zeros(Float64, Ms+1, Ms+1), #H
    zeros(ComplexF64, Ms+1, Ne), #psi
    zeros(ComplexF64, Ms+1, Ne), #phi
    zeros(Float64, Ms+1, Ms+1), #Γ
    MVector{C, Float64}(zeros(Float64, Ms+1)), #λ
    zeros(Float64, Ms+1, Ms+1), #dhdea
    zeros(Float64, Ms+1, Ms+1), #dhdv
    MVector{1, ComplexF64}(zero(ComplexF64)), #akl
    MVector{1, Float64}(0.0), #akk
    zeros(Float64, Ms+1,Ms+1), #DM
    zeros(Float64, Ne, Ms+1-Ne), #blk
    zeros(Float64, Ne, Ms+1-Ne), #blk2
    zeros(Float64, Ne*(Ms+1-Ne)), #Pb
    MVector{1, ComplexF64}(zero(ComplexF64)), #"phipsi"
    MVector{1, Float64}(0.0), #"Pbmaxest"
    zeros(Float64, N+2, 3),#vdot
    zeros(Float64, N+2, 3),#vtemp
    MVector{1, Float64}(0.0), #vscale
    MVector{1, Int64}(0), #exnum
    MVector{1, Int64}(0), #attnum
    MVector{1, Int64}(0), #nf
    zeros(ComplexF64, Ne, Ne), #temp_akl
    MVector{C, ComplexF64}(zeros(ComplexF64, C)), #uu
    zeros(Float64, Ms+1, Ne), #uuu
    zeros(Float64, Ms+1, Ms+1),   #temp_vm_gamma
    zeros(Float64, Ms+1, Ms+1), #Γ_hold
    MVector{C, Float64}(zeros(Float64, C)),#temp_prop_gamma
    zeros(ComplexF64, Ne, Ne), #temp_phi_psi
    zeros(Float64, tsteps), #storage_aop
    zeros(Float64, Ms+1, tsteps), #storage_op,
    zeros(Float64, 5, tsteps), #storage_e
    zeros(Float64, 6, tsteps), #storage_xno
    zeros(Float64, 6, tsteps), #storage_vno
    zeros(Float64, tsteps), #storage_temp
    zeros(Float64, tsteps), #storage_phop
    zeros(Float64, tsteps), #storage kinetic energy
    zeros(Float64, tsteps), #storage potential energy
    zeros(Float64, Int(round(tsteps/2))), #storage_state
    zeros(Float64, Int(round(tsteps/2))), #storage_deltaKNO
    zeros(Float64, Int(round(tsteps/2))), #storage_deltaKAu
    zeros(Int64, Int(round(tsteps/2))), #storage_hoptimes
    zeros(ComplexF64, (Ms+1)*Ne, tsteps), #storage_psi
    zeros(ComplexF64, (Ms+1)*Ne, tsteps), #storage_phi
    zeros(Float64, 396*3, tsteps), #storage_xau
    zeros(Float64, 396*3, tsteps)
    ) #storage_vau
elseif extensive_logging == false
    return Simulation{A, C, D}(
    zeros(Float64, N+2, 3), #x
    zeros(Float64, N+2, 3), #v
    MVector{3, Float64}(0.0, 0.0, 0.0),    #Δ_no
    zeros(Int64, N, 12), #nn_arr
    zeros(Float64, N+2, 3), #dhp_neutral
    zeros(Float64, N+2, 3), #dhp_ion
    zeros(Float64, N+2, 3), #dhp_coup
    MVector{3, Float64}(0.0, 0.0, 0.0), #hp
    zeros(Float64, N+2, 3), #F
    MVector{A, Int16}(zeros(Int16, Ne)), #surfp
    MVector{C, Int16}(zeros(Int16, Ms+1)), #occnum
    MVector{D, Int16}(zeros(Int16, Ms-Ne+1)), #surfh
    MVector{A, Int16}(zeros(Int16, Ne)), #surfpinit
    MVector{A, Int16}(zeros(Int16, Ne)), #surfpnew
    MVector{1, Float64}(10.0), # trajzmin
    MVector{1, Float64}(0.0), #trajtheta
    zeros(Float64, Ms+1, Ms+1), #H
    zeros(ComplexF64, Ms+1, Ne), #psi
    zeros(ComplexF64, Ms+1, Ne), #phi
    zeros(Float64, Ms+1, Ms+1), #Γ
    MVector{C, Float64}(zeros(Float64, Ms+1)), #λ
    zeros(Float64, Ms+1, Ms+1), #dhdea
    zeros(Float64, Ms+1, Ms+1), #dhdv
    MVector{1, ComplexF64}(zero(ComplexF64)), #akl
    MVector{1, Float64}(0.0), #akk
    zeros(Float64, Ms+1,Ms+1), #DM
    zeros(Float64, Ne, Ms+1-Ne), #blk
    zeros(Float64, Ne, Ms+1-Ne), #blk2
    zeros(Float64, Ne*(Ms+1-Ne)), #Pb
    MVector{1, ComplexF64}(zero(ComplexF64)), #"phipsi"
    MVector{1, Float64}(0.0), #"Pbmaxest"
    zeros(Float64, N+2, 3),#vdot
    zeros(Float64, N+2, 3),#vtemp
    MVector{1, Float64}(0.0), #vscale
    MVector{1, Int64}(0), #exnum
    MVector{1, Int64}(0), #attnum
    MVector{1, Int64}(0), #nf
    zeros(ComplexF64, Ne, Ne), #temp_akl
    MVector{C, ComplexF64}(zeros(ComplexF64, C)), #uu
    zeros(Float64, Ms+1, Ne), #uuu
    zeros(Float64, Ms+1, Ms+1),   #temp_vm_gamma
    zeros(Float64, Ms+1, Ms+1), #Γ_hold
    MVector{C, Float64}(zeros(Float64, C)),#temp_prop_gamma
    zeros(ComplexF64, Ne, Ne), #temp_phi_psi
    zeros(Float64, tsteps), #storage_aop
    zeros(Float64, Ms+1, tsteps), #storage_op,
    zeros(Float64, 5, tsteps), #storage_e
    zeros(Float64, 6, tsteps), #storage_xno
    zeros(Float64, 6, tsteps), #storage_vno
    zeros(Float64, tsteps), #storage_temp
    zeros(Float64, tsteps), #storage_phop
    zeros(Float64, tsteps), #storage kinetic energy
    zeros(Float64, tsteps), #storage potential energy
    zeros(Float64, Int(round(tsteps/2))), #storage_state
    zeros(Float64, Int(round(tsteps/2))), #storage_deltaKNO
    zeros(Float64, Int(round(tsteps/2))), #storage_deltaKAu
    zeros(Int64, Int(round(tsteps/2))) #storage_hoptimes
    )
end

@doc """
    simulation_constructor_x_v!(s)

    Initializes:
    s.x: position vector of N, O and Au atoms
    s.v: velocity vector of N, O and Au atoms
    s.Δ_no: Vector pointing from N to O.

"""
function simulation_constructor_x_v!(s::Simulation)

    vvib_no = sqrt(e_vib_i*2.0/(mass_arr[1]*mass_arr[2]/(mass_arr[1] + mass_arr[2])))# Initial vibrational velocity of NO
    vz_no = sqrt(e_trans_i*2.0/(mass_arr[1] + mass_arr[2]))# Initial z velocity of NO

    #initial position of NO
    x0no = zeros(Float64, 2, 3)

    xno_rand = 0.5#rand() #random number between (0, 1)
    xno = @. xno_rand * [cell[1], cell[2], 0.0]
    #xno[3] = 10.5 #Initial z position of NO = 10
    xno[3] = 1.5
    #Construct N and O position vector in simulation basis by transforming
    #from NO basis using angles theta and phi which are samples such that the
    #simulation basis unit vectors of NO fill out uniformly the 1-sphere
    #citation: http://corysimon.github.io/articles/uniformdistn-on-sphere/

    theta_no = 0.5#rand() # Initial orientation of NO, 0 = O-down, PI = N-down
    theta_no = acos(1.0 - 2.0 * theta_no)
    phi_no = 0.5#rand()
    phi_no = 2.0*π*phi_no

    x0no[1, :] = xno .+ r_no/2.0 * [sin(phi_no)*sin(theta_no), -cos(phi_no)*sin(theta_no), cos(theta_no)]
    x0no[2, :] = xno .- r_no/2.0 * [sin(phi_no)*sin(theta_no), -cos(phi_no)*sin(theta_no), cos(theta_no)]

    x0no[1, 3] = x0no[1, 3] + 0.2
    x0no[2, 2] = x0no[2, 2]  + 0.2
    xno = x0no[1,:] - x0no[2, :] #unit vector pointing from O to N
    xno = xno / norm(xno)

    #initial velocity of no
    v0no = zeros(Float64, 2, 3) #construct container for initial velocity of N and O
    phi_no = 0.0#2*π*rand() #direction of incidence, i.e. North, south east or west
    v0no[1, :] = vz_no * [sin(phi_no)*sin(phi_inc), -cos(phi_no)*sin(phi_inc), -cos(phi_inc)]
    v0no[2, :] = v0no[1, :]
    v0no[1, :] = v0no[1, :] + vvib_no * xno * mass_arr[2] / (mass_arr[1] + mass_arr[2])
    v0no[2, :] = v0no[2, :] - vvib_no * xno * mass_arr[1] / (mass_arr[1] + mass_arr[2])

    #add rotational velocity
    if e_rot_i != 0.0
        theta_no = 0.5#rand()
        theta_no = asin(2.0 * theta_no - 1.0) + π/2.0
        vrot_no = @. [-xno[3], 0.0, xno[1]] * cos(theta_no)/sqrt(xno[1]^2 + xno[3]^2)
        vrot_no = @. vrot_no + [xno[1], -(xno[1]^2 + xno[3]^2)/xno[2], xno[3]] *
                    sin(theta_no)/sqrt((xno[1]^2 + xno[3]^2) * (1 + (xno[1]^2 + xno[3]^2)/xno[2]^2))

        vvib_no = sqrt(2.0 * e_rot_i/(mass_arr[1] + mass_arr[2]))
        v0no[1, :] = @. v0no[1, :] + vvib_no * vrot_no * sqrt(mass_arr[2]/mass_arr[1])
        v0no[2, :] = @. v0no[2, :] - vvib_no * vrot_no * sqrt(mass_arr[1]/mass_arr[2])
    end
    #stack x0no and v0no to x_au0 and v_au0 to create x and v
    v0 = [-30.0*v0no; v_au0]
    #v0 = [v0no; v_au0]
    x0 = [x0no; x_au0]
    s.x .= x0
    s.v .= v0
    #construct difference vector between N and O for timestep t=0
    s.Δ_no .= x0[1, :] - x0[2, :]
end

@doc """
    simulation_constructor_nn!(s)

    Computes unit cell neighbors for each grid Atom
    Initializes: s.nn_arr
"""
function simulation_constructor_nn!(s::Simulation)
    nn_lattice!(s)
end

@doc """
    simulation_constructor_energy(s)

    Computes neutral, ionic and coupling energy for the NO diabatic Hamiltonian
    Initializes: s.hp (s.hp[1] ~ Neutral, s.hp[2] ~ Ionic, s.hp[3] ~ Coupling)
"""
function simulation_constructor_energy(s::Simulation)
    s.hp .= get_E_all(s)
end

@doc """
    simulation_constructor_force!(s)

    Computes partial derivatives of s.hp[i] with respect to the position vector of particle j
    Initializes: s.dhp_neutral, s.dhp_ionic, s.dhp_coup. Important: All partial derivatives
    have been multiplied by minus 1 to convert to forces.
"""
function simulation_constructor_force(s::Simulation)
    get_F_all(s)
end

@doc """
    simulation_constructor_x_300K!(s)

    Overwrites equilibrium position vectors of grid with position vectors of thermalized surface
"""
function simulation_constructor_x_300K!(s::Simulation)
    s.x[3:end, :] .= x_300K0
end
