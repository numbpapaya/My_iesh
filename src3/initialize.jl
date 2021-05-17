function initialize()
    #Initialize NO properties

    vvib_no = sqrt(e_vib_i*2.0/(mass_arr[1]*mass_arr[2]/(mass_arr[1] + mass_arr[2])))		# Initial vibrational velocity of NO
    vz_no = sqrt(e_trans_i*2.0/(mass_arr[1] + mass_arr[1]))				# Initial z velocity of NO

    #initial position of NO
    x0no = zeros(Float64, 2, 3)

    xno = rand() #temporary variable
    xno = @. xno * [cell[1], cell[2], 1.0]
    xno[3] = xno[3] + 10.5 # Initial z position of NO = 10

    theta_no = rand() # Initial orientation of NO, 0 = O-down, PI = N-down
    theta_no = asin(2.0 * theta_no - 1.0) + π/2.0
    phi_no = rand()
    phi_no = 2.0*π*phi_no



    x0no[1, :] = xno + r_no/2.0 * [sin(phi_no)*sin(theta_no), cos(phi_no)*sin(theta_no), cos(theta_no)]
    x0no[2, :] = xno - r_no/2.0 * [sin(phi_no)*sin(theta_no), cos(phi_no)*sin(theta_no), cos(theta_no)]

    xno = x0no[1,:] - x0no[2, :]
    xno = xno / sqrt(xno'*xno)

    #initial velocity of no
    v0no = zeros(Float64, 2, 3)
    phi_no = 0.0
    v0no[1, :] = vz_no * [sin(phi_no)*sin(phi_inc), -cos(phi_no)*sin(phi_inc), -cos(phi_inc)]
    v0no[2, :] = v0no[1, :]
    v0no[1, :] = v0no[1, :] + vvib_no * xno * mass_arr[2] / (mass_arr[1] + mass_arr[2])
    v0no[2, :] = v0no[2, :] - vvib_no * xno * mass_arr[2] / (mass_arr[1] + mass_arr[2])

    #add rotational velocity
    theta_no = rand()
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



    #initialize nearest neigbors

    return x0, v0
end
