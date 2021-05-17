#needs testing!!!!

function nn_loop_if!(r,  nn_arr,  i::Int,  j::Int)
    rb = U_sa * r
    s = @. Integer(round(rb/g))
    sint = Integer(s[1] * 9 + s[2] * 3 + s[3])
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
    return nn_arr
end

function nn_loop!(x::float_array,  nn_arr::Array{Int64, 2})
    for i in 1:N
        for j in 1:N
            temp = (x[j+2,  :] - x[i+2,  :])./cell
            temp = temp .- floor.(temp .+ 1.0/2)
            r = temp.*cell
            rmi = floor(sum(r.*r)*10000.0)/10000.0
            if rmi <= (g*1.01)^2
                nn_arr = nn_loop_if!(r,  nn_arr,  i,  j)
            end
        end
    end
end


function nn_lattice(x)
    nn_arr = zeros(Int64, N+2, 12)
    nn_loop!(x, nn_arr)
    return nn_arr
end

hmm= nn_lattice(s.x)
