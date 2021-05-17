export time_loop

function time_loop(n)
    global x, dr, F
    global v, storage_e
    global trajzmin, dhdea
    global trajtheta, dhdv, blk
    global eigvec_H, dhp
    global eigval_H, surfpinit
    global phi, psi, storage_aop, surfp, storage_op, storage_x, storage_v
    dr = x[1, :] - x[2, :]
    v_cond = (v[1, 3]*m_N + v[2, 3]*m_O)/(m_N + m_O)
    z_cond = (x[1, 3] + x[2, 3])/2.0
    if z_cond > z_end && v_cond > 0.0
        println("Boundary reached: End simulation")
        return
    end

    # Record minimial distance of z approach and orientational angle
    min_z_no = minimum(x[1:2, 3])
    if min_z_no <= trajzmin
        trajzmin = min_z_no
        trajtheta = acos((dr[3])/sqrt(dr' * dr))
    end

    # get state corresponding to current surface
    for i in 1:Ne
        phi[:, i] = eigvec_H[:, surfp[i]]
    end

    # Propagate electronic Hamiltonian
    psi = eigvec_H' * psi
    psi = @. exp(-1im * eigval_H * dt/hbar)

    psi = repeat(psi, Ne)
    psi = reshape(psi, (41, Ne))
    psi = @.psi * psi
    psi = eigvec_H * psi


    #recordings
    storage_aop[n] = sum(abs.(psi[1, :]).^2)
    storage_op_temp = sum(abs.((eigvec_H' * psi).^2), dims=2)
    storage_op[:, n] = storage_op_temp

    rtemp = sqrt(dr' * dr)
    storage_e[1, n] = 0.5 * (m_N + m_O)*sum(((m_N * v[1, :] + m_O * v[2, :])/(m_N + m_O)).^2)
    T_vib = 0.50 * μ * sum(((v[1, :] - v[2, :]).*(x[1, :] - x[2, :])/rtemp).^2)
    U_vib = F_n *(1.0 - exp(-δ_n * (rtemp -  r_0_NO)))^2
    E_vib = U_vib + T_vib
    storage_e[2, n] = E_vib
    T_tot = 0.5 * mass_arr[1] * sum(v[1, :].^2) + 0.5*mass_arr[2]*sum(v[2, :].^2)
    storage_e[3, 1] = T_tot - T_vib - storage_e[1, n]
    storage_e[4, 1] = sum(eigval_H[surfp]) - sum(eigval_H[surfpinit])

    #Propagate nuclear motion

    vdot = @. 1.0/m_spread * F
    x = @. x + v*dt + 0.5*vdot * dt^2
    vtemp = @. v + 0.5 *vdot * dt

    # Get adiabatic eigenvalues, eigenvectors
    hp = get_E_all(x) #Hamiltonian parameters
    dhp = get_F_all(x)

    # dhp = get_F_all(x, dhp)

    H = h0 + vm * hp[3] / sqrt(delta_E)
    H[1, 1] = H[1, 1] + hp[2] - hp[1]
    eig_vec_hold = eigvec_H

    eigval_H = eigvals(H)
    eigvec_H = eigvecs(H)

    for j in 1:Ms+1
        vHmaxloc = argmax(abs.(eig_vec_hold' * eigvec_H[:, j]))
        sign_b = sign(dot(eig_vec_hold[:, vHmaxloc], eigvec_H[:, j]))
        eigvec_H[:, j] = eigvec_H[:, j]*sign_b
    end
    # Get dH / d Ep matrices where Ep is parameter in Newns-Anderson Hamiltonian
    # (i.e. Ep = V or Ep = E_a = (EI-EN) )
    blk .= 0
    for j in 1:Ms+1
        dhdea[:, j] = @. eigvec_H[1, :] * eigvec_H[1, j]
    end
    dhdv = eigvec_H' * (vm * eigvec_H)/sqrt(delta_E)

    #calculate forces

    F .= 0
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
    v[399:end, :] .= 0.0
    F = -F

    v = @. vtemp + 0.5 / m_spread * F * dt

    #storage
    storage_x[1:3, n] = x[1, :]
    storage_x[4:6, n] = x[2, :]

    storage_v[1:3, n] = v[1, :]
    storage_v[4:6, n] = v[2, :]
    return
end
