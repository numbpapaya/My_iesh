
function simulation_init()
    A = Ne
    C = Ms+1
    D = Ms-Ne+1
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
    MVector{1, Float64}(10.0*Å), # trajzmin
    MVector{1, Float64}(0.0), #trajtheta
    zeros(Float64, Ms+1, Ms+1), #H
    zeros(ComplexF64, Ms+1, Ne), #psi
    zeros(ComplexF64, Ms+1, Ne), #phi
    zeros(Float64, Ms+1, Ms+1), #eigvec_H
    MVector{C, Float64}(zeros(Float64, Ms+1)), #eigval_H
    zeros(Float64, Ms+1, Ms+1), #dhdea
    zeros(Float64, Ms+1, Ms+1), #dhdv
    MVector{1, ComplexF64}(zero(ComplexF64)), #akl
    MVector{1, Float64}(0.0), #akk
    zeros(Float64, Ms+1,Ms+1), #DM
    zeros(Float64, Ne, Ms+1-Ne), #blk
    zeros(Float64, Ms+1, Ms+1), #blk2
    zeros(Float64, Ne*(Ms+1-Ne)), #Pb
    zeros(Float64, tsteps), #storage_aop
    zeros(Float64, Ms+1, tsteps), #storage_op,
    zeros(Float64, 5, tsteps), #storage_e
    zeros(Float64, 6, tsteps), #storage_xno
    zeros(Float64, 6, tsteps), #storage_vno
    zeros(Float64, tsteps), #storage_temp
    zeros(Float64, tsteps), #storage_phop
    MVector{1, ComplexF64}(zero(ComplexF64)), #"phipsi"
    MVector{1, Float64}(0.0), #"Pbmaxest"
    zeros(Float64, N+2, 3),#vdot
    zeros(Float64, N+2, 3),#vtemp
    MVector{1, Float64}(0.0), #vscale
    zeros(Float64, tsteps), #KEt
    zeros(Float64, tsteps), #PEt
    zeros(Float64, Int(round(tsteps/2))),
    zeros(Float64, Int(round(tsteps/2))),
    zeros(Float64, Int(round(tsteps/2))),
    zeros(Int64, Int(round(tsteps/2))),
    MVector{1, Int64}(0),
    MVector{1, Int64}(0),
    MVector{1, Int64}(0),
    zeros(ComplexF64, (Ms+1)*Ne, tsteps), #storage_psi
    zeros(ComplexF64, (Ms+1)*Ne, tsteps), #storage_phi
    zeros(Float64, 396*6, tsteps), #storage_xau
    zeros(Float64, 396*6, tsteps) #storage_vau
    )
end


function simulation_constructor_x_v!(s::Simulation)

    vvib_no = sqrt(e_vib_i*2.0/(mass_arr[1]*mass_arr[2]/(mass_arr[1] + mass_arr[2])))		# Initial vibrational velocity of NO
    vz_no = sqrt(e_trans_i*2.0/(mass_arr[1] + mass_arr[2]))				# Initial z velocity of NO

    #initial position of NO
    x0no = zeros(Float64, 2, 3)

    xno = rand() #temporary variable
    #xno = 0.1#rand() #temporary variable
    xno = @. xno * [cell[1], cell[2], 1.0*Å]
    #xno[3] = xno[3] + 4.5*Å # Initial z position of NO = 10
    xno[3] = xno[3] + 10.5*Å
    theta_no = rand() # Initial orientation of NO, 0 = O-down, PI = N-down
    #theta_no = 0.1#rand() # Initial orientation of NO, 0 = O-down, PI = N-down
    theta_no = asin(2.0 * theta_no - 1.0) + π/2.0
    phi_no = rand()
    #phi_no = 0.1#rand()
    phi_no = 2.0*π*phi_no



    x0no[1, :] = xno .+ r_no/2.0 * [sin(phi_no)*sin(theta_no), cos(phi_no)*sin(theta_no), cos(theta_no)]
    x0no[2, :] = xno .- r_no/2.0 * [sin(phi_no)*sin(theta_no), cos(phi_no)*sin(theta_no), cos(theta_no)]
    #x0no[1, 3] = x0no[1, 3] +1Å
    xno = x0no[1,:] - x0no[2, :]
    xno = xno / norm(xno)

    #initial velocity of no
    v0no = zeros(Float64, 2, 3)
    phi_no = 0.0
    v0no[1, :] = vz_no * [sin(phi_no)*sin(phi_inc), -cos(phi_no)*sin(phi_inc), -cos(phi_inc)]
    v0no[2, :] = v0no[1, :]
    v0no[1, :] = v0no[1, :] + vvib_no * xno * mass_arr[2] / (mass_arr[1] + mass_arr[2])
    v0no[2, :] = v0no[2, :] - vvib_no * xno * mass_arr[1] / (mass_arr[1] + mass_arr[2])

    #add rotational velocity
    theta_no = rand()
    #theta_no = 0.1#rand()
    theta_no = asin(2.0 * theta_no - 1.0) + π/2.0
    vrot_no = @. [-xno[3], 0.0, xno[1]] * cos(theta_no)/sqrt(xno[1]^2 + xno[3]^2)
    vrot_no = @. vrot_no + [xno[1], -(xno[1]^2 + xno[3]^2)/xno[2], xno[3]] *
                sin(theta_no)/sqrt((xno[1]^2 + xno[3]^2) * (1 + (xno[1]^2 + xno[3]^2)/xno[2]^2))

    vvib_no = sqrt(2.0 * e_rot_i/(mass_arr[1] + mass_arr[2]))
    v0no[1, :] = @. v0no[1, :] + vvib_no * vrot_no * sqrt(mass_arr[2]/mass_arr[1])
    v0no[2, :] = @. v0no[2, :] - vvib_no * vrot_no * sqrt(mass_arr[1]/mass_arr[2])
    #stack information to x_au0 and v_au0
    v0 = [v0no; v_au0]
    x0 = [x0no; x_au0]
    s.x .= x0
    s.v .= v0
    s.Δ_no .= x0[1, :] - x0[2, :]
end

function simulation_constructor_nn!(s::Simulation)
    nn_lattice!(s)
end


function simulation_constructor_energy(s::Simulation)
    s.hp .= get_E_all(s)
end

function simulation_constructor_force(s::Simulation)
    get_F_all(s)
end
