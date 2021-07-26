export nn_lattice!

@inline @inbounds function nn_loop_if!(r,  nn_arr,  i::Integer,  j::Integer, rb, ss)
    mul!(rb, U_sa, r)
    ss .= Integer(round(rb./g))
    sint = Integer(ss[1] * 9 + ss[2] * 3 + ss[3]) #nice hack, it works, but why? figure & table  1, Begbie 1947
    if sint == 4
        nn_arr[i,  1] = j
    elseif sint == 2
        nn_arr[i,  2] = j
    elseif sint == 10
        nn_arr[i,  3] = j
    elseif sint == -8
        nn_arr[i,  4] = j
    elseif sint == 12
        nn_arr[i,  5] = j
    elseif sint == 6
        nn_arr[i,  6] = j
    elseif sint == -4
        nn_arr[i,  7] = j
    elseif sint == -2
        nn_arr[i,  8] = j
    elseif sint == -10
        nn_arr[i,  9] = j
    elseif sint == 8
        nn_arr[i,  10] = j
    elseif sint == -12
        nn_arr[i,  11] = j
    elseif sint == -6
        nn_arr[i,  12] = j
    end
end

@inline @inbounds function nn_loop!(x,  nn_arr, temp_vec, r, rmi, rb, s)
    for i in 1:N
        for j in 1:N
            @views temp_vec .= (x[j+2,  :] .- x[i+2,  :])./cell
            temp_vec .= temp_vec .- floor.(temp_vec .+ 1.0/2.0)
            r .= temp_vec .* cell
            r_dot = norm(r)
            if r_dot <= (g*1.01)
                nn_loop_if!(r,  nn_arr,   i,  j, rb, s)
            end
        end
    end
end


@inline function nn_lattice!(s::Simulation)
    temp_vec = zeros(Float64, 3)
    r = zeros(Float64, 3)
    rmi = 0.0
    rb = zeros(Float64, 3)
    ss = zeros(Float64, 3)
    nn_arr = nn_loop!(s.x, s.nn_arr, temp_vec, r, rmi, rb, ss)

end
