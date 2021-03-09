#needs testing!!!!

function nn(x):
    k = 13
    nn_data = similar(x)
    nn_data = transpose(x)
    kdtree = KDTree(nn_data)
    idxs, dists = knn(kdtree, nn_data, k)
    idxs_temp = zeros(Int8, k, Integer(N))
    dists_temp = zeros(Float64, k, Integer(N))
    idxs_temp = hcat(idxs...)
    dists_temp = hcat(dists...)
    return idxs_temp, dists_temp
end
