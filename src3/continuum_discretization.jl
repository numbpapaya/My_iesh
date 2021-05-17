function burkey_cantrell()
    e_diabat = zeros(Float64, Ms)
    h0 = zeros(Float64, Ms+1, Ms+1)
    vm = zeros(Float64, Ms+1, Ms+1)
    gauss = zeros(Float64, Int(Ms/2), Int(Ms/2))
    for j in 1:(Int(Ms/2)-1)
        gauss[j, j+1] = Float64(j)/sqrt(4.0*Float64(j)^2 - 1.0)
        gauss[j+1, j] = gauss[j, j+1]
    end

    eigvals_gauss = eigvals(gauss)
    eigvecs_gauss = eigvecs(gauss)

    for j in 2:(Int(Ms/2) + 1)
        h0[j, j] = delta_E/4.0*eigvals_gauss[j-1] - delta_E/4.0
        h0[j+Int(Ms/2), j+Int(Ms/2)] = delta_E/4.0*eigvals_gauss[j-1] + delta_E/4.0
        e_diabat[j-1] = h0[j, j]
        e_diabat[j + Int(Ms/2) - 1] = h0[j+Int(Ms/2), j+Int(Ms/2)]
    end

    vm[1, 2:Int(Ms/2) + 1] = @. sqrt(delta_E/4.0 * (2.0*eigvecs_gauss[1, :]^2))
    vm[2:Int(Ms/2) + 1, 1] = @. sqrt(delta_E/4.0 * (2.0*eigvecs_gauss[1, :]^2))
    vm[1, (2 + Int(Ms/2)):(Ms + 1)] = @. sqrt(delta_E/4.0 * (2.0*eigvecs_gauss[1, :]^2))
    vm[(2 + Int(Ms/2)):(Ms + 1), 1] = @. sqrt(delta_E/4.0 * (2.0*eigvecs_gauss[1, :]^2))

    return h0, e_diabat, vm
end
