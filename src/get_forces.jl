export get_F_all

function get_F_all(s::Simulation)
    F_au = zeros(Float64, N+2, 3)
    F_au_n_n= zeros(Float64, N+2, 3)
    F_au_o_n= zeros(Float64, N+2, 3)
    F_n_o_n_temp = MMatrix{2, 3, Float64, 6}(zeros(2, 3))
    F_n_o_n= zeros(Float64, N+2, 3)
    F_neutral= zeros(Float64, N+2, 3)

    F_au_n_i = zeros(Float64, N+2, 3)
    F_au_o_i= zeros(Float64, N+2, 3)
    F_n_o_i_temp = MMatrix{2, 3, Float64, 6}(zeros(2, 3))
    F_n_o_i = zeros(Float64, N+2, 3)
    F_image_i_temp = MMatrix{2, 3, Float64, 6}(zeros(2, 3))
    F_image_i = zeros(Float64, N+2, 3)
    F_ion = zeros(Float64, N+2, 3)


    F_au_n_c = zeros(Float64, N+2, 3)
    F_au_o_c = zeros(Float64, N+2, 3)
    F_coupling = zeros(Float64, N+2, 3)
    get_F_all_inner!(s, F_au, F_au_n_n , F_au_o_n,F_n_o_n_temp ,F_n_o_n,F_neutral,F_au_n_i ,F_au_o_i,F_n_o_i_temp ,F_n_o_i ,F_image_i_temp ,F_image_i ,F_ion ,F_au_n_c ,F_au_o_c ,F_coupling)
end


function get_F_all_inner!(s::Simulation, F_au, F_au_n_n , F_au_o_n,F_n_o_n_temp ,F_n_o_n,F_neutral,F_au_n_i ,F_au_o_i,F_n_o_i_temp ,F_n_o_i ,F_image_i_temp ,F_image_i ,F_ion ,F_au_n_c ,F_au_o_c ,F_coupling
)

    F_au .= get_F_au_au(s)
    F_au_n_n .= get_F_au_n_neutral(s)
    F_au_o_n .= get_F_au_o_neutral(s)
    F_n_o_n_temp .= get_F_n_o_neutral(s)
    F_n_o_n .= [F_n_o_n_temp; zeros(Float64, N, 3)]
    F_neutral .= @. F_au_n_n + F_au_o_n + F_n_o_n - F_au #check #add half force of gold-gold?

    F_au_n_i .= get_F_au_n_ion(s) #check
    F_au_o_i.= get_F_au_o_ion(s) #check
    F_n_o_i_temp .= get_F_n_o_ion(s)
    F_n_o_i .= [F_n_o_i_temp; zeros(Float64, N, 3)]
    F_image_i_temp .= get_F_image_ion(s)
    F_image_i .= [F_image_i_temp; zeros(Float64, N, 3)]
    F_ion .= @. F_au_n_i + F_au_o_i + F_n_o_i + F_image_i - F_au #check #add half force of gold-gold?


    F_au_n_c .= get_F_au_n_coup(s)
    F_au_o_c .= get_F_au_o_coup(s)
    F_coupling .= @. F_au_n_c + F_au_o_c

    #diabatic matrix elements of NO-Au Hamiltonian
    s.dhp_neutral .= -F_neutral
    s.dhp_ion .= -F_ion
    s.dhp_coup .= -F_coupling
end
#---------------------------Au-N Coupling----------------------------------
@inline function get_F_au_n_coup_loop!(x::float_array, F::float_array, i,
     delta_xn, cutoff_test, bool_cutoff_test, s::Simulation)::float_array
    @views delta_xn =  x[i+2, :] - x[1, :]
    delta_xn[1] = @.delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = @.delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test .= abs.(delta_xn) .- cutoff
    bool_cutoff_test .= cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xn = norm(delta_xn)
        if norm_delta_xn <= cutoff
            term1 = exp(coup_β_N*norm_delta_xn)
            term2 = (1.0 + coup_b_N*term1)^2 * norm_delta_xn
            term3 = coup_a_N * coup_b_N * coup_β_N/term2

            dV_dri = @. -term1*delta_xn*term3
            dV_drn = -dV_dri
            @views F[1, :] .= F[1, :] .- dV_drn
            @views F[i+2, :] .= F[i+2, :] .- dV_dri
        end
    end
    return F
end

@inbounds @inline function get_F_au_n_coup(s::Simulation)::float_array
    F = zeros(Float64, N+2, 3)
    delta_xn = MVector{3, Float64}(zeros(Float64, 3))
    cutoff_test = MVector{3, Float64}(zeros(Float64, 3))
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        F = get_F_au_n_coup_loop!(s.x, F, i, delta_xn, cutoff_test, bool_cutoff_test, s)
    end
    return F
end
#---------------------------Au-O Coupling----------------------------------
@inline function get_F_au_o_coup_loop!(x::float_array, F::float_array, i, delta_xo,
    cutoff_test, bool_cutoff_test, s::Simulation)::float_array
    @views delta_xo = x[i+2, :] .- x[2, :]
    delta_xo[1] = delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test .= abs.(delta_xo) .- cutoff
    bool_cutoff_test .= cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xo = norm(delta_xo)
        if norm_delta_xo <= cutoff
            term1 = exp(coup_β_O*norm_delta_xo)
            term2 = (1.0 + coup_b_O*term1)^2 * norm_delta_xo
            term3 = coup_a_O * coup_b_O * coup_β_O/term2

            dV_dri = -term1*delta_xo*term3
            dV_dro = -dV_dri
            @views F[2, :] .= F[2, :] .- dV_dro
            @views F[i+2, :] .= F[i+2, :] .- dV_dri
        end
    end
    return F
end

@inline @inbounds function get_F_au_o_coup(s::Simulation)::float_array
    F = zeros(Float64, N+2, 3)
    delta_xo = MVector{3, Float64}(zeros(Float64, 3))
    cutoff_test = MVector{3, Float64}(zeros(Float64, 3))
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        F = get_F_au_o_coup_loop!(s.x, F, i, delta_xo, cutoff_test, bool_cutoff_test, s)
    end
    return F
end

#---------------------------Ionic N O ----------------------------------

function get_F_n_o_ion(s::Simulation)::MMatrix{2, 3, Float64, 6}
    F = MMatrix{2, 3, Float64, 6}(zeros(2, 3))
    norm_r_no = norm(s.Δ_no)
    term = 2*F_i *δ_i * (exp(-2.0 * δ_i * (norm_r_no - rNO_e_i)) - exp(-δ_i*(norm_r_no - rNO_e_i)))
    F[1, :] .=  term/norm_r_no .* r_no
    F[2, :] .= - F[1, :]
    return F
end
#---------------------------Ionic Au O ----------------------------------
@inline function get_F_au_o_ion_loop(x::Array{Float64,2}, F::float_array, i::Int,
     delta_xo, cutoff_test, bool_cutoff_test, s::Simulation)::float_array
    @views delta_xo .=  x[i+2, :] .- x[2, :]
    delta_xo[1] = delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test .= abs.(delta_xo) .- cutoff
    bool_cutoff_test .= cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xo = norm(delta_xo)
        if norm_delta_xo <= cutoff
            term = -A_n * exp(-α_n*norm_delta_xo) * α_n/norm_delta_xo
            dV_dri = term * delta_xo
            dV_dro = -dV_dri
            F[2, :] .= view(F, 2, :) .- dV_dro
            F[i+2, :] .= view(F, i+2, :) .- dV_dri
        end
    end
    return F
end

@inline @inbounds function get_F_au_o_ion(s::Simulation)::float_array #check
    F = zeros(Float64, N+2, 3)
    delta_xo = MVector{3, Float64}(zeros(Float64, 3))
    cutoff_test = MVector{3, Float64}(zeros(Float64, 3))
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        F = get_F_au_o_ion_loop(s.x, F, i, delta_xo, cutoff_test, bool_cutoff_test, s)
    end
    return F
end

#---------------------------Ionic Au N ----------------------------------
@inline @inbounds function get_F_au_n_ion_loop(x::Array{Float64,2}, F::float_array,
     i::Int, delta_xn, cutoff_test, bool_cutoff_test, s::Simulation)::float_array #check
    @views delta_xn .=  x[i+2, :] .- x[1, :]
    delta_xn[1] = delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])
    cutoff_test .=  abs.(delta_xn) .- cutoff
    bool_cutoff_test .=  cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xn = norm(delta_xn)
        if norm_delta_xn <= cutoff
            norm_r_no = norm(s.Δ_no)
            r_diff = norm_delta_xn - rN_sur_e_i
            term1 = s.x[2, 3] - s.x[1, 3]
            term2 = exp(-β_i*r_diff)
            term3 = exp(-2.0 * β_i * r_diff)
            term4 = norm_r_no^4
            #println(term2)
            dV_xN_1 = 4.0*term2*(s.Δ_no[1])*term1^2
            dV_xN_1 = dV_xN_1/term4
            dV_xN_4 = 2.0*term3*delta_xn[1]*β_i
            dV_xN_4 = dV_xN_4/norm_delta_xn
            dV_xN_5 = 2.0*term2*delta_xn[1]*term1*term1*β_i
            dV_xN_5 = dV_xN_5/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_xN = B_i*(dV_xN_1+dV_xN_4-dV_xN_5)

            dV_yN_1 = 4.0*term2*(s.Δ_no[2])*term1*term1
            dV_yN_1 = dV_yN_1/term4
            dV_yN_4 = 2.0*term3*delta_xn[2]*β_i
            dV_yN_4 = dV_yN_4/norm_delta_xn
            dV_yN_5 = 2.0*term2*delta_xn[2]*term1*term1*β_i
            dV_yN_5 = dV_yN_5/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_yN = B_i*(dV_yN_1+dV_yN_4-dV_yN_5)

            dV_zN_1 = 4.0*term2*(s.Δ_no[3])*term1*term1
            dV_zN_1 = dV_zN_1/term4
            dV_zN_2 = 4.0*term2*term1
            dV_zN_2 = dV_zN_2/(norm_r_no*norm_r_no)
            dV_zN_4 = 2.0*term3*delta_xn[3]*β_i
            dV_zN_4 = dV_zN_4/norm_delta_xn
            dV_zN_5 = 2.0*term2*delta_xn[3]*term1*term1*β_i
            dV_zN_5 = dV_zN_5/(norm_delta_xn*norm_r_no*norm_r_no)
            dV_zN = B_i*(dV_zN_1+dV_zN_2+dV_zN_4-dV_zN_5)

            dV_xO_1 = 4.0*term2*(s.Δ_no[1])*term1*term1
            dV_xO_1 = dV_xO_1/term4
            dV_xO = B_i*(-dV_xO_1)

            dV_yO_1 = 4.0*term2*(s.Δ_no[2])*term1*term1
            dV_yO_1 = dV_yO_1/term4
            dV_yO = B_i*(-dV_yO_1)

            dV_zO_1 = 4.0*term2*(s.Δ_no[3])*term1*term1
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

            dVc_xN_1 = 4.0*term2*(s.Δ_no[1])*term1*term1
            dVc_xN_1 = dVc_xN_1/term4
            dVc_xN = B_i*(dVc_xN_1)

            dVc_yN_1 = 4.0*term2*(s.Δ_no[2])*term1*term1
            dVc_yN_1 = dVc_yN_1/term4
            dVc_yN = B_i*(dVc_yN_1)

            dVc_zN_1 = 4.0*term2*(s.x[1, 3] - s.x[2, 3])*term1*term1
            dVc_zN_1 = dVc_zN_1/term4
            dVc_zN_2 = 4.0*term2*term1
            dVc_zN_2 = dVc_zN_2/(norm_r_no*norm_r_no)
            dVc_zN = B_i*(dVc_zN_1+dVc_zN_2)

            dVc_xO_1 = 4.0*term2*(s.Δ_no[1])*term1*term1
            dVc_xO_1 = dVc_xO_1/term4
            dVc_xO = B_i*(-dVc_xO_1)

            dVc_yO_1 = 4.0*term2*(s.Δ_no[2])*term1*term1
            dVc_yO_1 = dVc_yO_1/term4
            dVc_yO = B_i*(-dVc_yO_1)

            dVc_zO_1 = 4.0*term2*(s.Δ_no[3])*term1*term1
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

@inline function get_F_au_n_ion(s::Simulation)::float_array
    F = zeros(Float64, N+2, 3)
    delta_xn = MVector{3, Float64}(zeros(Float64, 3))
    cutoff_test = MVector{3, Float64}(zeros(Float64, 3))
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        F = get_F_au_n_ion_loop(s.x, F, i, delta_xn, cutoff_test, bool_cutoff_test, s)
    end
    return F
end



#---------------------------Ionic Image Potential----------------------------------
function get_F_image_ion(s::Simulation)::MMatrix{2, 3, Float64, 6}  #check
    F = MMatrix{2, 3, Float64, 6}(zeros(2, 3))
    @views zcom = zcom_NO(s.x[1:2, :])
    term = (D_i*(zcom - z0_i))/(C_i^2 + (zcom-z0_i)^2)^(1.5)
    F[1, 3] = term * mass_arr[1]/(mass_arr[1] + mass_arr[2])
    F[2, 3] = term * mass_arr[2]/(mass_arr[1] + mass_arr[2])
    return -F
end
#---------------------------Ionic N-O----------------------------------
function get_F_n_o_ion(s::Simulation)::MMatrix{2, 3, Float64, 6}
    F = MMatrix{2, 3, Float64, 6}(undef)
    norm_r_no = norm(s.Δ_no)
    term = 2*F_i*δ_i*exp(-δ_i*(norm_r_no - rNO_e_i))*(1-exp(-δ_i*(norm_r_no - rNO_e_i)))
    F[1, :] .= term/norm_r_no .* s.Δ_no
    F[2, :] .= -F[1, :]
    return -F
end

#---------------------------Neutral N-O----------------------------------

@inline @inbounds function get_F_n_o_neutral(s::Simulation)::MMatrix{2, 3, Float64, 6}
    F = MMatrix{2, 3, Float64, 6}(undef)
    norm_r_no = norm(s.Δ_no)
    rdiff = norm_r_no - r_0_NO
    term = 2*F_n*δ_n*(exp(-2.0*δ_n * rdiff) - exp(-δ_n * rdiff))
    #F[1, 1] = -term * (s.x[2, 1] - s.x[1, 1]) / norm_r_no
    F[1, :] .= term/norm_r_no .* s.Δ_no
    F[2, :] .= -F[1, :]
    return F
end
#---------------------------Neutral AU_O----------------------------------
@inline function get_F_au_o_neutral_loop(x::float_array, i::Int, F::float_array,
     delta_xo, cutoff_test, bool_cutoff_test, s::Simulation)::float_array
    @views delta_xo .= x[i+2, :] .- x[2, :]
    delta_xo[1] = delta_xo[1] - cell[1]*round(delta_xo[1]/cell[1])
    delta_xo[2] = delta_xo[2] - cell[2]*round(delta_xo[2]/cell[2])

    cutoff_test .= abs.(delta_xo) .- cutoff
    bool_cutoff_test .= cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xo = norm(delta_xo)
        if norm_delta_xo <= cutoff
            term = -A_n * exp(-α_n * norm_delta_xo) * α_n/norm_delta_xo
            dV_dxi = term * delta_xo
            dV_dxo = -dV_dxi

            F[2,:] .= view(F, 2, :) .- dV_dxo
            F[i+2,:] .= view(F, i+2, :) .- dV_dxi
        end
    end
    return F
end

@inbounds @inline function get_F_au_o_neutral(s::Simulation)::float_array
    F = zeros(Float64, N+2, 3)
    delta_xo = MVector{3, Float64}(zeros(Float64, 3))
    cutoff_test = MVector{3, Float64}(zeros(Float64, 3))
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        F = get_F_au_o_neutral_loop(s.x, i, F, delta_xo, cutoff_test, bool_cutoff_test, s)
    end
    return F
end

#---------------------------Neutral AU_N----------------------------------

@inline function get_F_au_n_neutral_loop(x::Array{Float64,2}, i::Int, F::float_array,
     delta_xn, cutoff_test, bool_cutoff_test, s::Simulation)::float_array #maybe check
    @views delta_xn .=  x[i+2, :] .- x[1, :]
    delta_xn[1] = delta_xn[1] - cell[1]*round(delta_xn[1]/cell[1])
    delta_xn[2] = delta_xn[2] - cell[2]*round(delta_xn[2]/cell[2])

    cutoff_test .= abs.(delta_xn) .- cutoff
    bool_cutoff_test .=  cutoff_test .<= 0.0
    sum_bool = sum(bool_cutoff_test)
    if sum_bool != zero(sum_bool)
        norm_delta_xn = norm(delta_xn)
        if norm_delta_xn <= cutoff
            term = -B_n * exp(-β_n * norm_delta_xn) * β_n/norm_delta_xn
            dV_dxi =  term * delta_xn
            dV_dxN = -dV_dxi

            F[1,:] .= view(F, 1, :) .- dV_dxN
            F[i+2,:] .= view(F, i+2, :) .- dV_dxi
        end
    end
    return F
end

@inbounds @inline function get_F_au_n_neutral(s::Simulation)::float_array #check
    F = zeros(Float64, N+2, 3)
    delta_xn = MVector{3, Float64}(zeros(Float64, 3))
    cutoff_test = MVector{3, Float64}(zeros(Float64, 3))
    bool_cutoff_test = BitArray{1}(undef, 3)
    for i in 1:N
        F = get_F_au_n_neutral_loop(s.x, i, F, delta_xn, cutoff_test, bool_cutoff_test, s)
    end
    return F
end

#---------------------------GOLD AU AU LATTICE----------------------------------
@inline function F_au_au_loop_if(x::float_array, F::float_array, i::Int, j::Int,
     xm, xi, temp,r, qf_1, s::Simulation)::float_array
    if s.nn_arr[i, j] != 0
        m = s.nn_arr[i, j] + 2
        @views xm .= x[m, :]
        @views xi .= x[i+2, :]
        temp .= xm .- xi
        temp .= temp./cell
        temp .= temp .- floor.(temp .+ 0.5)
        @views r .=  temp .* cell .- r0[:, j]
        if j == 1 || j== 7
            mul!(qf_1, d1_new, r)
            F[i+2, :] .= view(F, i+2, :) .- qf_1
        elseif j == 2 || j== 8
            mul!(qf_1, d2_new, r)
            F[i+2, :] .= view(F, i+2, :) .- qf_1
        elseif j == 3 || j== 9
            mul!(qf_1, d3_new, r)
            F[i+2, :] .= view(F, i+2, :) .- qf_1
        elseif j == 4 || j== 10
            mul!(qf_1, d4_new, r)
            F[i+2, :] .= view(F, i+2, :) .- qf_1
        elseif j == 5|| j== 11
            mul!(qf_1, d5_new, r)
            F[i+2, :] .= view(F, i+2, :) .- qf_1
        elseif j == 6 || j== 12
            mul!(qf_1, d6_new, r)
            F[i+2, :] .= view(F, i+2, :) .- qf_1
        end
    end
    return F
end
@inline @inbounds function F_au_au_loop(x::float_array, F::float_array,
     s::Simulation)::float_array #check
    xm = MVector{3, Float64}(zeros(Float64, 3))
    xi = MVector{3, Float64}(zeros(Float64, 3))
    temp = MVector{3, Float64}(zeros(Float64, 3))
    r = MVector{3, Float64}(zeros(Float64, 3))
    qf_1 = MVector{3, Float64}(zeros(Float64, 3))
    for i in 1:N
        for j in 1:12
            F .= F_au_au_loop_if(x, F, i, j, xm, xi, temp,r, qf_1, s)
        end
    end
    return F
end

@inline function get_F_au_au(s::Simulation)::float_array
    F = zeros(Float64, N+2, 3)
    F = F_au_au_loop(s.x, F, s)
    return F
end
