export get_F_all

function get_F_all(x::Array{Float64,2})
    dhp = zeros(Float64, 3, 3*(N+2))
    #Units are STILL not clear
    F = zeros(Float64, N, 3)
    get_f = get_F_au_au(x, F)
    F_au_au = [zeros(Float64, 2, 3); get_f]

    F_au_n_n = get_F_au_n_neutral(x)
    F_au_o_n = get_F_au_o_neutral(x)
    F_n_o_n = [get_F_n_o_neutral(x[1:2,:]); zeros(Float64, N, 3)]
    F_neutral = @. F_au_n_n + F_au_o_n + F_n_o_n + F_au_au #add half force of gold-gold?

    F_au_n_i = get_F_au_n_ion(x)
    F_au_o_i= get_F_au_o_ion(x)
    F_n_o_i = [get_F_n_o_ion(x[1:2,:]); zeros(Float64, N, 3)]
    zcom = zcom_NO(x[1:2, :])
    F_image_i = [get_F_image_ion(zcom); zeros(Float64, N, 3)]
    F_ion = @. F_au_n_i + F_au_o_i + F_n_o_i + F_image_i + F_au_au  #add half force of gold-gold?


    F_au_n_c = get_F_au_n_coup(x)
    F_au_o_c = get_F_au_o_coup(x)
    F_coupling = @. F_au_n_c + F_au_o_c

    #diabatic matrix elements of NO-Au Hamiltonian
    dhp[1, :] = vec(F_neutral')
    dhp[2, :] = vec(F_ion')
    dhp[3, :] = vec(F_coupling')
    return -dhp
end
#---------------------------Au-O Coupling----------------------------------
function get_F_au_n_coup_loop!(x::Array{Float64,2}, F::Array{Float64,2}, i::Int)
    delta_xn = @. x[i+2, :] - x[1, :]
    delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test = @. abs(delta_xn) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xn = sqrt(delta_xn' * delta_xn)
        if norm_delta_xn <= cutoff
            term1 = exp(coup_β_N*norm_delta_xn)
            term2 = (1.0 + coup_b_N*term1)^2 * norm_delta_xn
            term3 = coup_a_N * coup_b_N * coup_β_N/term2

            dV_dri = @. -term1*delta_xn*term3
            dV_drn = -dV_dri
            F[1, :] = F[1, :] - dV_drn
            F[i+2, :] = F[i+2, :] - dV_dri
        end
    end
    return F
end

function get_F_au_n_coup(x::Array{Float64,2})
    F = zeros(Float64, N+2, 3)
    for i in 1:N
        F = get_F_au_n_coup_loop!(x, F, i)
    end
    return F
end
#---------------------------Au-O Coupling----------------------------------
function get_F_au_o_coup_loop!(x::Array{Float64,2}, F::Array{Float64,2}, i::Int)
    delta_xo =  x[i+2, :] .- x[2, :]
    delta_xo[1] = @.delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = @.delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test = @. abs(delta_xo) - cutoff
    bool_cutoff_test =  cutoff_test .<= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xo = sqrt(delta_xo' * delta_xo)
        if norm_delta_xo <= cutoff
            term1 = exp(coup_β_O*norm_delta_xo)
            term2 = (1.0 + coup_b_O*term1)^2 * norm_delta_xo
            term3 = coup_a_O * coup_b_O * coup_β_O/term2

            dV_dri = @. -term1*delta_xo*term3
            dV_dro = -dV_dri
            F[2, :] = F[2, :] - dV_dro
            F[i+2, :] = F[i+2, :] - dV_dri
        end
    end
    return F
end

function get_F_au_o_coup(x::Array{Float64,2})
    F = zeros(Float64, N+2, 3)
    for i in 1:N
        F = get_F_au_o_coup_loop!(x, F, i)
    end
    return F
end

#---------------------------Ionic N O ----------------------------------

function get_F_n_o_ion(x_no::Array{Float64,2})
    F = zeros(Float64, 2, 3)
    r_no = x_no[1, :] - x_no[2,:]
    norm_r_no = sqrt(r_no' * r_no)
    term = 2*F_i *δ_i * (exp(-2.0 * δ_i * (norm_r_no - rNO_e_i)) - exp(-δ_i*(norm_r_no - rNO_e_i)))
    F[1, :] = @. term * r_no/norm_r_no
    F[2, :] = - F[1, :]
    return F
end

#---------------------------Ionic Au O ----------------------------------
function get_F_au_o_ion_loop!(x::Array{Float64,2}, F::Array{Float64,2}, i::Int)
    delta_xo =  x[i+2, :] .- x[2, :]
    delta_xo[1] = @.delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = @.delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test = @. abs(delta_xo) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xo = sqrt(delta_xo' * delta_xo)
        if norm_delta_xo <= cutoff
            term = -A_n * exp(-α_n*norm_delta_xo) * α_n/norm_delta_xo
            dV_dri = term * delta_xo
            dV_dro = -dV_dri
            F[2, :] = @. F[2, :] - dV_dro
            F[i+2, :] = @. F[i+2, :] - dV_dri
        end
    end
    return F
end

function get_F_au_o_ion(x::Array{Float64,2})
    F = zeros(Float64, N+2, 3)
    for i in 1:N
        F = get_F_au_o_ion_loop!(x, F, i)
    end
    return F
end



#---------------------------Ionic Au N ----------------------------------

function get_F_au_n_ion_loop!(x::Array{Float64,2}, F::Array{Float64,2}, i::Int)
    delta_xn =  x[i+2, :] .- x[1, :]
    delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test = @. abs(delta_xn) - cutoff
    bool_cutoff_test =  cutoff_test .<= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xn = sqrt(delta_xn' * delta_xn)
        if norm_delta_xn <= cutoff
            r_no = x[1, :] - x[2,:]
            norm_r_no = sqrt(r_no' * r_no)
            r_diff = norm_delta_xn - rN_sur_e_i
            term1 = r_no[3]
            term2 = exp(-β_i*r_diff)
            term3 = exp(-2.0 * β_i * r_diff)
            term4 = norm_r_no^4

            dV_xN_1 = 4.0*term2*(x[1, 1]-x[2, 1])*term1*term1
            dV_xN_1 = dV_xN_1/term4
            dV_xN_4 = 2.0*term3*delta_xn[1]*β_i
            dV_xN_4 = dV_xN_4/norm_delta_xn
            dV_xN_5 = 2.0*term2*delta_xn[1]*term1*term1*β_i
            dV_xN_5 = dV_xN_5/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_xN = B_i*(dV_xN_1+dV_xN_4-dV_xN_5)

            dV_yN_1 = 4.0*term2*(x[1, 2]-x[2, 2])*term1*term1
            dV_yN_1 = dV_yN_1/term4
            dV_yN_4 = 2.0*term3*delta_xn[2]*β_i
            dV_yN_4 = dV_yN_4/norm_delta_xn
            dV_yN_5 = 2.0*term2*delta_xn[2]*term1*term1*β_i
            dV_yN_5 = dV_yN_5/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_yN = B_i*(dV_yN_1+dV_yN_4-dV_yN_5)

            dV_zN_1 = 4.0*term2*(x[1, 3]-x[2, 3])*term1*term1
            dV_zN_1 = dV_zN_1/term4
            dV_zN_2 = 4.0*term2*term1
            dV_zN_2 = dV_zN_2/(norm_r_no*norm_r_no)
            dV_zN_4 = 2.0*term3*delta_xn[3]*β_i
            dV_zN_4 = dV_zN_4/norm_delta_xn
            dV_zN_5 = 2.0*term2*delta_xn[3]*term1*term1*β_i
            dV_zN_5 = dV_zN_5/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_zN = B_i*(dV_zN_1+dV_zN_2+dV_zN_4-dV_zN_5)

            dV_xO_1 = 4.0*term2*(x[1, 1]-x[1, 2])*term1*term1
            dV_xO_1 = dV_xO_1/term4
            dV_xO = B_i*(-dV_xO_1)

            dV_yO_1 = 4.0*term2*(x[1, 1]-x[2, 2])*term1*term1
            dV_yO_1 = dV_yO_1/term4
            dV_yO = B_i*(-dV_yO_1)

            dV_zO_1 = 4.0*term2*(x[1, 3]-x[2, 3])*term1*term1
            dV_zO_1 = dV_zO_1/term4
            dV_zO_2 = 4.0*term2*term1
            dV_zO_2 = dV_zO_2/(norm_r_no*norm_r_no)
            dV_zO = B_i*(-dV_zO_1-dV_zO_2)

            dV_xi_3 = 2.0*term3*delta_xn[1]*β_i
            dV_xi_3 = dV_xi_3/norm_delta_xn
            dV_xi_4 = 2.0*term2*delta_xn[1]*term1*term1*β_i
            dV_xi_4 = dV_xi_4/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_xi = B_i*(-dV_xi_3+dV_xi_4)

            dV_yi_3 = 2.0*term3*delta_xn[2]*β_i
            dV_yi_3 = dV_yi_3/norm_delta_xn
            dV_yi_4 = 2.0*term2*delta_xn[2]*term1*term1*β_i
            dV_yi_4 = dV_yi_4/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_yi = B_i*(-dV_yi_3+dV_yi_4)

            dV_zi_3 = 2.0*term3*delta_xn[3]*β_i
            dV_zi_3 = dV_zi_3/norm_delta_xn
            dV_zi_4 = 2.0*term2*delta_xn[3]*term1*term1*β_i
            dV_zi_4 = dV_zi_4/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_zi = B_i*(-dV_zi_3+dV_zi_4)

            #------   Calculating derivatives of cutoff function

            r_diff = cutoff - rN_sur_e_i

            term2 = exp(-β_i*r_diff)
            term3 = exp(-2.0*β_i*r_diff)

            dVc_xN_1 = 4.0*term2*(x[1, 1]-x[2, 1])*term1*term1
            dVc_xN_1 = dVc_xN_1/term4
            dVc_xN = B_i*(dVc_xN_1)

            dVc_yN_1 = 4.0*term2*(x[1, 2]-x[2, 2])*term1*term1
            dVc_yN_1 = dVc_yN_1/term4
            dVc_yN = B_i*(dVc_yN_1)

            dVc_zN_1 = 4.0*term2*(x[1, 3]-x[2, 3])*term1*term1
            dVc_zN_1 = dVc_zN_1/term4
            dVc_zN_2 = 4.0*term2*term1
            dVc_zN_2 = dVc_zN_2/(norm_r_no*norm_r_no)
            dVc_zN = B_i*(dVc_zN_1+dVc_zN_2)

            dVc_xO_1 = 4.0*term2*(x[1, 1]-x[2, 1])*term1*term1
            dVc_xO_1 = dVc_xO_1/term4
            dVc_xO = B_i*(-dVc_xO_1)

            dVc_yO_1 = 4.0*term2*(x[1, 2]-x[2, 2])*term1*term1
            dVc_yO_1 = dVc_yO_1/term4
            dVc_yO = B_i*(-dVc_yO_1)

            dVc_zO_1 = 4.0*term2*(x[1, 3]-x[2, 3])*term1*term1
            dVc_zO_1 = dVc_zO_1/term4
            dVc_zO_2 = 4.0*term2*term1
            dVc_zO_2 = dVc_zO_2/(norm_r_no*norm_r_no)
            dVc_zO = B_i*(-dVc_zO_1-dVc_zO_2)

            dVc_xi = 0.0
            dVc_yi = 0.0
            dVc_zi = 0.0

            F[1, 1] = F[1, 1] - dV_xN + dVc_xN
            F[1, 2] = F[1, 2] - dV_yN + dVc_yN
            F[1, 3] = F[1, 3] - dV_zN + dVc_zN

            F[2, 1] = F[2, 1] - dV_xO + dVc_xO
            F[2, 2] = F[2, 2] - dV_yO + dVc_yO
            F[2, 3] = F[2, 3] - dV_zO + dVc_zO

            F[i+2, 1] = F[i+2, 1] - dV_xi + dVc_xi
            F[i+2, 2] = F[i+2, 2] - dV_yi + dVc_yi
            F[i+2, 3] = F[i+2, 3] - dV_zi + dVc_zi
        end
    end
    return F
end

function get_F_au_n_ion(x::Array{Float64,2})
    F = zeros(Float64, N+2, 3)
    for i in 1:N
        F = get_F_au_n_ion_loop!(x, F, i)
    end
    return F
end

#---------------------------Ionic Image Potential----------------------------------
function get_F_image_ion(zcom::Float64)
    F = zeros(Float64, 2, 3)
    term = (D_i*(zcom - z0_i))/(C_i^2 + (zcom-z0_i)^2)^(1.5)
    F[1, 3] = term * mass_arr[1]/(mass_arr[1] + mass_arr[2])
    F[2, 3] = term * mass_arr[2]/(mass_arr[1] + mass_arr[2])
    return -F
end
#---------------------------Ionic N-O----------------------------------
function get_F_n_o_ion(x_no::Array{Float64,2})
    F = zeros(Float64, 2, 3)
    r_no = x_no[1, :] - x_no[2,:]
    norm_r_no = sqrt(r_no' * r_no)
    term = 2*F_i*δ_i*exp(-δ_i*(norm_r_no - rNO_e_i))*(1-exp(-δ_i*(norm_r_no - rNO_e_i)))
    F[1, :] = @. term * r_no/norm_r_no
    F[2, :] = -F[1, :]
    return -F
end

#---------------------------Neutral N-O----------------------------------

function get_F_n_o_neutral(x_no::Array{Float64,2})
    F = zeros(Float64, 2, 3)
    r_no = x_no[1, :] - x_no[2,:]
    norm_r_no = sqrt(r_no' * r_no)
    term = 2*F_n*δ_n*exp(-δ_n*(norm_r_no - r_0_NO))*(1-exp(-δ_n*(norm_r_no - r_0_NO)))
    F[1, :] = @. term * r_no/norm_r_no
    F[2, :] = -F[1, :]
    return -F
end

#---------------------------Neutral AU_O----------------------------------
function get_F_au_o_neutral_loop!(x_all::Array{Float64,2}, i::Int, F::Array{Float64,2})
    delta_xo =  x_all[i+2, :] .- x_all[2, :]
    delta_xo[1] = @. delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = @. delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test = @. abs(delta_xo) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xo = sqrt(delta_xo' * delta_xo)
        if norm_delta_xo <= cutoff
            term = -A_n * exp(-α_n * norm_delta_xo) * α_n/norm_delta_xo
            dV_dxi = @. term * delta_xo
            dV_dxO = -dV_dxi

            F[2,:] = F[2,:] - dV_dxO
            F[i+2,:] = F[i+2,:] - dV_dxi
        end
    end
    return F
end

function get_F_au_o_neutral(x_all::Array{Float64,2})
    F::Array{Float64,2} = zeros(Float64, N+2, 3) #Force acting on N and all Au atoms
    for i in 1:N
        F = get_F_au_o_neutral_loop!(x_all, i, F)
    end
    return F
end


#---------------------------Neutral AU_N----------------------------------
function get_F_au_n_neutral_loop!(x_all::Array{Float64,2}, i::Int, F::Array{Float64,2})
    delta_xn =  x_all[i+2, :] .- x_all[1, :]
    delta_xn[1] = @. delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @. delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test = @. abs(delta_xn) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xn = sqrt(delta_xn' * delta_xn)
        if norm_delta_xn <= cutoff
            term = -B_n * exp(-β_n * norm_delta_xn) * β_n/norm_delta_xn
            dV_dxi = @. term * delta_xn
            dV_dxN = -dV_dxi

            F[1,:] = F[1,:] + dV_dxN
            F[i+2,:] = F[i+2,:] + dV_dxi
        end
    end
    return F
end

function get_F_au_n_neutral(x_all::Array{Float64,2})
    F::Array{Float64,2} = zeros(Float64, N+2, 3) #Force acting on N and all Au atoms
    for i in 1:N
        F = get_F_au_n_neutral_loop!(x_all, i, F)
    end
    return -F
end



#---------------------------GOLD AU AU LATTICE----------------------------------
function F_au_au_loop_if!(x::Array{Float64,2}, F::Array{Float64,2}, i::Int, j::Int)
    if nn[i, j] != 0
        m::Int = nn[i, j]+2
        xm = x[m, :]
        xi = x[i+2, :]
        temp = xm .- xi
        temp .= temp./cell
        temp_floor = floor.(temp .+ 0.5)
        temp .=  temp .- temp_floor
        r = temp .* cell .- r0[:, j]

        if j == 1 || j== 7
            F[i,:] = F[i,:]  .+ (d1_new + d1_new') * r
        elseif j == 2 || j== 8
            F[i,:] = F[i,:]  .+  (d2_new + d2_new') * r
        elseif j == 3 || j== 9
            F[i,:] = F[i,:]  .+  (d3_new + d3_new') * r
        elseif j == 4 || j== 10
            F[i,:] = F[i,:]  .+  (d4_new + d4_new') * r
        elseif j == 5|| j== 11
            F[i,:] = F[i,:]  .+  (d5_new + d5_new') * r
        elseif j == 6 || j== 12
            F[i,:] = F[i,:]  .+  (d6_new + d6_new') * r
        end
    end
    return F
end

function F_au_au_loop!(x_au::Array{Float64,2}, F::Array{Float64,2})
    for i in 1:N
        for j in 1:12
            F = F_au_au_loop_if!(x_au, F, i, j)
        end
    end
    return F
end

function get_F_au_au(x_au::Array{Float64,2}, F::Array{Float64,2})
    #Force acting on all N AU atoms
    F = F_au_au_loop!(x_au, F)
    return -F
end
