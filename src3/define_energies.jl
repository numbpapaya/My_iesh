export get_E_all
#ENERGIES
function get_E_all(x)
    #Units are STILL not clear
    V = zero(Float64)
    E_au_au = get_V_au_au(x, V)
    E_coup = get_E_coup(x)
    E_neutral = get_E_neutral(x) + E_au_au
    E_ion = get_E_ion(x) + E_au_au
    hp = [E_neutral, E_ion, E_coup]
    return hp
end

#---------------------------Coupling----------------------------------

function get_E_coup(x::Array{Float64,2})
    E_au_n_coup = get_E_au_n_coup(x)
    E_au_o_coup = get_E_au_o_coup(x)
    return E_au_n_coup + E_au_o_coup
end

#---------------------------Au-N Coupling----------------------------------
function get_E_au_n_coup_loop!(x, V, i)
    delta_xn = @. x[i+2, :] - x[1, :]
    delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test = @. abs(delta_xn) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xn = sqrt(delta_xn' * delta_xn)
        if norm_delta_xn <= cutoff
            V_iterate = coup_a_N * (
            1.0/ (1.0 + coup_b_N * exp(coup_β_N*norm_delta_xn)) - coup_cutoff_N
            )
            V = V + V_iterate
        end
    end
    return V
end

function get_E_au_n_coup(x)
    V = 0.0
    for i in 1:N
        V = get_E_au_n_coup_loop!(x, V, i)
    end
    return V
end
#---------------------------Au-O Coupling----------------------------------
function get_E_au_o_coup_loop!(x, V, i)
    delta_xo = @. x[i+2, :] - x[2, :]
    delta_xo[1] = @.delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = @.delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test = @. abs(delta_xo) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xo = sqrt(delta_xo' * delta_xo)
        if norm_delta_xo <= cutoff
            V_coup = coup_a_O /(1.0 + coup_b_O * exp(coup_β_O*norm_delta_xo))
            V_cut = coup_cutoff_O
            V = V + V_coup - V_cut
        end
    end
    return V
end

function get_E_au_o_coup(x)
    V = 0.0
    for i in 1:N
        V = get_E_au_o_coup_loop!(x, V, i)
    end
    return V
end
#---------------------------Ionic----------------------------------

function get_E_ion(x)
    E_n_o_ion= get_E_n_o_ion(x[1:2, :])
    E_au_n_ion = get_E_au_n_ion(x)
    E_au_o_ion = get_E_au_o_ion(x)
    zcom = zcom_NO(x[1:2, :])
    E_image = get_E_image_ionic(zcom)
    return E_n_o_ion + E_au_n_ion+ E_au_o_ion + E_image + KI
end

#---------------------------Ionic N O ----------------------------------

function get_E_n_o_ion(x_no::Array{Float64,2})
    r_no = x_no[1, :] - x_no[2,:]
    norm_r_no = sqrt(r_no' * r_no)
    return F_i * (1 - exp(-δ_i*(norm_r_no - rNO_e_i)))^2
end

#---------------------------Ionic Au O ----------------------------------
function get_E_au_o_ion_loop!(x::Array{Float64,2}, V::Float64, i::Int)
    delta_xo = @. x[i+2, :] - x[2, :]
    delta_xo[1] = @.delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = @.delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test = @. abs(delta_xo) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xo = sqrt(delta_xo' * delta_xo)
        if norm_delta_xo <= cutoff
            V_iterate = A_n * (
            exp(-α_i * norm_delta_xo)
            - exp_alpha_i_cutoff
            )
            V = V + V_iterate
        end
    end
    return V
end

function get_E_au_o_ion(x::Array{Float64,2})
    V = 0.0
    for i in 1:N
        V = get_E_au_o_ion_loop!(x, V, i)
    end
    return V
end


#---------------------------Ionic Au N ----------------------------------

function get_E_au_n_ion_loop!(x::Array{Float64,2}, V::Float64, i::Int)
    delta_xn = @. x[i+2, :] - x[1, :]
    delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test = @. abs(delta_xn) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xn = sqrt(delta_xn' * delta_xn)
        if norm_delta_xn <= cutoff
            r_no = x[1, :] - x[2,:]
            norm_r_no = sqrt(r_no' * r_no)
            cos_theta = (x[1, 3] - x[2, 3])/norm_r_no
            V_iterate = @. B_i * (
            exp(-2*β_i *(norm_delta_xn - rN_sur_e_i)) -
            exp_beta_i_cutoff1 -
            2*cos_theta^2 * exp(-β_i *(norm_delta_xn - rN_sur_e_i))
            - exp_beta_i_cutoff2
            )
            V = V + V_iterate
        end
    end
    return V
end

function get_E_au_n_ion(x::Array{Float64,2})
    V = 0.0
    for i in 1:N
        V = get_E_au_n_ion_loop!(x, V, i)
    end
    return V
end

#---------------------------Ionic Image Potential----------------------------------
function get_E_image_ionic(zcom::Float64)
    return -D_i/sqrt(C_i^2 + ((zcom - z0_i))^2)
end
#---------------------------Ionic N-O----------------------------------
function get_E_n_o_ion(x_no::Array{Float64,2})
    r_no = x_no[1, :] - x_no[2,:]
    norm_r_no = sqrt(r_no' * r_no)
    return F_i * (1 - exp(-δ_i*(norm_r_no - rNO_e_i)))^2
end

#---------------------------Neutral----------------------------------

function get_E_neutral(x::Array{Float64,2})
    E_n_o_neutral = get_E_n_o_neutral(x[1:2, :])
    E_au_n_neutral = get_E_au_n_neutral(x)
    E_au_o_neutral = get_E_au_o_neutral(x)
    return E_n_o_neutral + E_au_n_neutral + E_au_o_neutral
end


#---------------------------Neutral N-O----------------------------------

function get_E_n_o_neutral(x_no::Array{Float64,2})
    r_no = @. x_no[1, :] - x_no[2,:]
    norm_r_no = sqrt(r_no' * r_no)
    return F_n * (1 - exp(-δ_n*(norm_r_no - r_0_NO)))^2
end



#---------------------------Neutral AU_O----------------------------------
function get_E_au_o_neutral_loop!(x, i, V)
    delta_xn = @. x[i+2, :] - x[2, :]
    delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test = @. abs(delta_xn) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xn = sqrt(delta_xn' * delta_xn)
        if norm_delta_xn <= cutoff
            V = V + A_n * (exp(-α_n*norm_delta_xn) - exp_alpha_n_cutoff)
        end
    end
    return V
end

function get_E_au_o_neutral(x)
    V = 0.0
    for i in 1:N
        V = get_E_au_o_neutral_loop!(x, i, V)
    end
    return V
end

#---------------------------Neutral AU_N----------------------------------
# function get_E_au_n_neutral_loop!(x::Array{Float64,2}, i::Int, V::Float64)
#     delta_xn = @. x[i+2, :] - x[1, :]
#     delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
#     delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])
#
#     cutoff_test = @. abs(delta_xn) - cutoff
#     bool_cutoff_test = @. cutoff_test <= 0.0
#     if sum(bool_cutoff_test) != 0.0
#         norm_delta_xn = sqrt(delta_xn' * delta_xn)
#         if norm_delta_xn <= cutoff
#             V = V + B_n * (exp(-β_n*norm_delta_xn) - exp_beta_n_cutoff)
#         end
#     end
#     return V
# end
#
# function get_E_au_n_neutral(x::Array{Float64,2})
#     V = 0.0
#     for i in 1:N
#         V = get_E_au_n_neutral_loop!(x, i, V)
#     end
#     return V
# end
function get_E_au_n_neutral_loop!(x::Array{Float64,2}, i::Int, V::Float64)
    delta_xn = @. x[i+2, :] - x[1, :]
    delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test = @. abs(delta_xn) - cutoff
    bool_cutoff_test = @. cutoff_test <= 0.0
    if sum(bool_cutoff_test) != 0.0
        norm_delta_xn = sqrt(delta_xn' * delta_xn)
        if norm_delta_xn <= cutoff
            V = V + B_n * (exp(-β_n*norm_delta_xn) - exp_beta_n_cutoff)
        end
    end
    return V
end

function get_E_au_n_neutral(x::Array{Float64,2})
    V = 0.0
    for i in 1:N
        V = get_E_au_n_neutral_loop!(x, i, V)
    end
    return V
end

#---------------------------GOLD AU AU LATTICE----------------------------------
function V_au_au_loop_if!(x::Array{Float64,2}, V::Float64, i::Int, j::Int)
    if nn[i, j] != 0
        m::Int = nn[i, j]+2
        xm = x[m, :]
        xi = x[i+2, :]
        temp = xm .- xi
        temp .= temp./cell
#        @views temp =  (x[m, :] .- x[i+2, :])./cell
        temp_floor = floor.(temp .+ 0.5)
        temp =  temp .- temp_floor
        r =  temp .* cell .- r0[:, j]

        if j == 1 || j== 7
            V = V + dot(r, d1_new * r)
        elseif j == 2 || j== 8
            V = V + dot(r, d2_new * r)
        elseif j == 3 || j== 9
            V = V + dot(r, d3_new * r)
        elseif j == 4 || j== 10
            V = V + dot(r, d4_new * r)
        elseif j == 5|| j== 11
            V = V + dot(r, d5_new * r)
        elseif j == 6 || j== 12
            V = V + dot(r, d6_new * r)
        end
    end
    return V
end
function V_au_au_loop!(x::Array{Float64,2}, V::Float64)
    for i in 1:N
        for j in 1:12
            V = V_au_au_loop_if!(x, V, i, j)
        end
    end
    return V
end

function get_V_au_au(x::Array{Float64,2}, V)
    V = V_au_au_loop!(x, V)
    return V/4.0
end
