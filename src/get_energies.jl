
#ENERGIES
@inline function get_E_all(s::Simulation)::MVector{3, Float64}
    E_au_au = get_V_au_au(s)
    E_coup = get_E_coup(s)
    E_neutral = get_E_neutral(s) + E_au_au
    E_ion = get_E_ion(s) + E_au_au
    hp = MVector{3, Float64}(E_neutral, E_ion, E_coup)
    return hp
end

#---------------------------Coupling----------------------------------

function get_E_coup(s::Simulation)::Float64 #check
    E_au_n_coup = get_E_au_n_coup(s)
    E_au_o_coup = get_E_au_o_coup(s)
    return E_au_n_coup + E_au_o_coup
end

#---------------------------Au-N Coupling----------------------------------
@inline function get_E_au_n_coup_loop!(x, V, i, delta_xn, cutoff_test,
    bool_cutoff_test, s::Simulation)::Float64
    @views delta_xn =  x[i+2, :] - x[1, :]
    delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test .= abs.(delta_xn) .- cutoff
    bool_cutoff_test .= cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xn = norm(delta_xn)
        if norm_delta_xn <= cutoff
            V_iterate = coup_a_N / (1.0 + coup_b_N * exp(coup_β_N*norm_delta_xn))
            V_cut = coup_cutoff_N
            V = V + V_iterate - V_cut
        end
    end
    return V
end

@inbounds @inline function get_E_au_n_coup(s::Simulation)::Float64
    V = 0.0
    delta_xn = Vector{Float64}(undef, 3)
    cutoff_test = Vector{Float64}(undef, 3)
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        V = get_E_au_n_coup_loop!(s.x, V, i, delta_xn, cutoff_test, bool_cutoff_test, s)
    end
    return V
end
#---------------------------Au-O Coupling----------------------------------
@inline function get_E_au_o_coup_loop!(x, V, i, delta_xo, cutoff_test,
     bool_cutoff_test, s::Simulation)::Float64
    @views delta_xo = x[i+2, :] .- x[2, :]
    delta_xo[1] = delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test .= abs.(delta_xo) .- cutoff
    bool_cutoff_test .= cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xo = norm(delta_xo)
        if norm_delta_xo <= cutoff
            V_coup = coup_a_O /(1.0 + coup_b_O * exp(coup_β_O*norm_delta_xo))
            V_cut = coup_cutoff_O
            V = V + V_coup - V_cut
        end
    end
    return V
end

@inline @inbounds function get_E_au_o_coup(s::Simulation)::Float64
    V = 0.0
    delta_xo = Vector{Float64}(undef, 3)
    cutoff_test = Vector{Float64}(undef, 3)
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        V = get_E_au_o_coup_loop!(s.x, V, i, delta_xo, cutoff_test, bool_cutoff_test, s)
    end
    return V
end
#---------------------------Ionic----------------------------------

function get_E_ion(s::Simulation)::Float64  #check
    E_n_o_ion= get_E_n_o_ion(s)
    E_au_n_ion = get_E_au_n_ion(s)
    E_au_o_ion = get_E_au_o_ion(s)
    E_image = get_E_image_ionic(s)
    return E_n_o_ion + E_au_n_ion+ E_au_o_ion + E_image + KI
end

#---------------------------Ionic N O ----------------------------------

function get_E_n_o_ion(s::Simulation)::Float64
    norm_r_no = norm(s.Δ_no)
    return F_i * (1 - exp(-δ_i*(norm_r_no - rNO_e_i)))^2
end

#---------------------------Ionic Au O ----------------------------------
@inline function get_E_au_o_ion_loop(x::Array{Float64,2}, V::Float64, i::Int,
    delta_xo, cutoff_test, bool_cutoff_test, s::Simulation)::Float64
    @views delta_xo .=  x[i+2, :] .- x[2, :]
    delta_xo[1] = delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test .= abs.(delta_xo) .- cutoff
    bool_cutoff_test .= cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xo = norm(delta_xo)
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

@inline @inbounds function get_E_au_o_ion(s::Simulation)::Float64
    V = 0.0
    delta_xo = Vector{Float64}(undef, 3)
    cutoff_test = Vector{Float64}(undef, 3)
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        V = get_E_au_o_ion_loop(s.x, V, i, delta_xo, cutoff_test, bool_cutoff_test, s)
    end
    return V
end

#---------------------------Ionic Au N ----------------------------------
@inline @inbounds function get_E_au_n_ion_loop(x::Array{Float64,2}, V::Float64,
     i::Int, delta_xn, cutoff_test, bool_cutoff_test, s::Simulation)::Float64
    @views delta_xn .=  x[i+2, :] .- x[1, :]
    delta_xn[1] = delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test .=  abs.(delta_xn) .- cutoff
    bool_cutoff_test .=  cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xn = norm(delta_xn)
        if norm_delta_xn <= cutoff
            norm_Δ_no = norm(s.Δ_no)
            cos_theta = (s.Δ_no[3])/norm_Δ_no
            V_iterate =  B_i * (
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

@inline function get_E_au_n_ion(s::Simulation)::Float64
    V = 0.0
    delta_xn = Vector{Float64}(undef, 3)
    cutoff_test = Vector{Float64}(undef, 3)
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        V = get_E_au_n_ion_loop(s.x, V, i, delta_xn, cutoff_test, bool_cutoff_test, s)
    end
    return V
end

#---------------------------Ionic Image Potential----------------------------------
@inline function get_E_image_ionic(s::Simulation)::Float64
    @views zcom = zcom_NO(s.x[1:2, :])
    return -D_i/sqrt(C_i^2 + ((zcom - z0_i))^2)
end
#---------------------------Ionic N-O----------------------------------
@inline function get_E_n_o_ion(s::Simulation)::Float64
    @views r_no = s.x[1, :] - s.x[2,:]
    norm_r_no = norm(r_no)
    return F_i * (1 - exp(-δ_i*(norm_r_no - rNO_e_i)))^2
end
#---------------------------Neutral----------------------------------

@inline function get_E_neutral(s::Simulation)::Float64  #check
    E_n_o_neutral = get_E_n_o_neutral(s) #check
    E_au_n_neutral = get_E_au_n_neutral(s) #check
    E_au_o_neutral = get_E_au_o_neutral(s) #check
    return E_n_o_neutral + E_au_n_neutral + E_au_o_neutral
end


#---------------------------Neutral N-O----------------------------------

@inline function get_E_n_o_neutral(s::Simulation)::Float64 #its good
    norm_r_no = norm(s.Δ_no)
    return F_n * (1 - exp(-δ_n*(norm_r_no - r_0_NO)))^2
end



#---------------------------Neutral AU_O----------------------------------
@inline function get_E_au_o_neutral_loop!(x, i, V, delta_xn, cutoff_test,
     bool_cutoff_test, s::Simulation)::Float64 #i think its ok, small numerical instability
    @views delta_xn .= x[i+2, :] .- x[2, :]
    delta_xn[1] = delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test .= abs.(delta_xn) .- cutoff
    bool_cutoff_test .= cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xn = norm(delta_xn)
        if norm_delta_xn <= cutoff
            V = V + A_n * (exp(-α_n*norm_delta_xn) - exp_alpha_n_cutoff)
        end
    end
    return V
end

@inbounds @inline function get_E_au_o_neutral(s::Simulation)::Float64
    V = 0.0
    delta_xn = zeros(Float64, 3)
    cutoff_test = zeros(Float64, 3)
    bool_cutoff_test = zeros(Bool, 3)
    for i in 1:N
        V = get_E_au_o_neutral_loop!(s.x, i, V, delta_xn, cutoff_test, bool_cutoff_test, s)
    end
    return V
end
#---------------------------Neutral AU_N----------------------------------

@inline function get_E_au_n_neutral_loop(x::Array{Float64,2}, i::Int, V::Float64,
    delta_xn, cutoff_test, bool_cutoff_test, s::Simulation)::Float64
    @views delta_xn .=  x[i+2, :] .- x[1, :]
    delta_xn[1] = delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test .= abs.(delta_xn) .- cutoff
    bool_cutoff_test .=  cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    V_temp = 0.0
    if sum_bool != zero(sum_bool)
        norm_delta_xn = norm(delta_xn)
        if norm_delta_xn <= cutoff
            V_temp = B_n * (exp(-β_n*norm_delta_xn) - exp_beta_n_cutoff)
        end
    end
    return V_temp
end

@inbounds @inline function get_E_au_n_neutral(s::Simulation)::Float64
    V = 0.0
    delta_xn = zeros(Float64, 3)
    cutoff_test = zeros(Float64, 3)
    bool_cutoff_test = zeros(Bool, 3)
    for i in 1:N
        V_temp =  get_E_au_n_neutral_loop(s.x, i, V, delta_xn, cutoff_test, bool_cutoff_test, s)
        V = V + V_temp
    end
    return V
end

#---------------------------GOLD AU AU LATTICE----------------------------------
@inline function V_au_au_loop_if(x::float_array, V::Float64, i::Int, j::Int,
     xm, xi, temp,r, qf_1, qf_2, s::Simulation)::Float64

    if s.nn_arr[i, j] != 0
        m = s.nn_arr[i, j] + 2
        @views xm .= x[m, :]
        @views xi .= x[i+2, :]
        temp .= xm .- xi
        temp .= temp./cell
        temp .= temp .- floor.(temp .+ 0.5)
        @views r .=  temp .* cell .- r0[:, j]
        if j == 1 || j== 7
            mul!(qf_2, d1_new, r)
            qf_1 = dot(r, qf_2)
            V = V + qf_1
        elseif j == 2 || j== 8
            mul!(qf_2, d2_new, r)
            qf_1 = dot(r, qf_2)
            V = V + qf_1
        elseif j == 3 || j== 9
            mul!(qf_2, d3_new, r)
            qf_1 = dot(r, qf_2)
            V = V + qf_1
        elseif j == 4 || j== 10
            mul!(qf_2, d4_new, r)
            qf_1 = dot(r, qf_2)
            V = V + qf_1
        elseif j == 5|| j== 11
            mul!(qf_2, d5_new, r)
            qf_1 = dot(r, qf_2)
            V = V + qf_1
        elseif j == 6 || j== 12
            mul!(qf_2, d6_new, r)
            qf_1 = dot(r, qf_2)
            V = V + qf_1
        end
    end
    return V
end
@inline @inbounds function V_au_au_loop(x::float_array, V::Float64, s::Simulation)::Float64
    xm = zeros(Float64, 3)
    xi = zeros(Float64, 3)
    temp = zeros(Float64, 3)
    r = zeros(Float64, 3)
    qf_1 = 0.0
    qf_2 = zeros(Float64, 3)
    for i in 1:N
        for j in 1:12
            V = V_au_au_loop_if(x, V, i, j, xm, xi, temp,r, qf_1, qf_2, s)
        end
    end
    return V
end

@inline function get_V_au_au(s::Simulation)::Float64    #check
    V = 0.0
    V = V_au_au_loop(s.x, V, s)
    return V*0.25 #why is here a factor of 1/4 instead of 1/2?
end
