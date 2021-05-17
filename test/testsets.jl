using My_iesh
using Test


# @testset "My_iesh.jl" begin
#     # Write your tests here.
#     @test f(2, 1) == 7
#     @test f(2, 3) == 13
# end
#
# @testset "oh hallo" begin
#     @test derivative_of_f(2, 1) == 2
end
#test NearestNeighbors
# plt3d= Plots.plot(testo3[:,1],testo3[:,2], testo3[:,3],
#      seriestype=:scatter, markersize = 7)
# display(plt3d)
#
#
# testo2 = similar(testo)
# testo2[:, 1] = testo[:, 1] .- testo[1,1]
# testo2[:, 2] = testo[:, 2] .- testo[1,2]
# testo2[:, 3] = testo[:, 3] .- testo[1,3]
#
# testo4 = similar(testo)
# for i in 1:13
#     testo3[i, :] = U' * testo2[i, :]
# end
# #
# # x = test_u2[:, 1]
# #
# # y = test_u2[:, 2]
# #
# # z = test_u2[:, 3]
# #
# # plt3d= Plots.plot(x,y, z,
# #      seriestype=:scatter, markersize = 7)
# # display(plt3d)
# #
# # theta = -π/4
# # rx = [[1 0 0]; [0 cos(theta) sin(theta)]; [0 -sin(theta) cos(theta)]]
# # ry = [[cos(theta) 0 -sin(theta)]; [0 1 0]; [sin(theta) 0 cos(theta)]]
# # r111 = rx * ry
# #
# #
# miu = zeros(Int32, 12)
# U = permutedims([[-1 0 1]/sqrt(2);[1 -2 1]/sqrt(6);[-1 -1 -1]/sqrt(3)])
# # for i in 1:13
#     for j in 2:13
#         temp = (testo[j, :] - testo[i, :])./abc
#         temp = temp .- floor.(temp .+ 1.0/2)
#         r = temp.*abc
#         rmi = floor(sum(r.*r)*10000)/10000
#         if rmi <= g^2
#             rb = U * r
#             # s = Integer.(sign.(rb) .* round.(abs.(rb)/g))
#             s = @. Integer(round(rb/g))
#             sint = Integer(s[1] * 9 + s[2] * 3 + s[3])
#             println(sint)
#             if sint == 4
#                 miu[1] = j
#             elseif sint == 2
#                 miu[2] = j
#             elseif sint == 10
#                 miu[3] = j
#             elseif sint == -8
#                 miu[4] = j
#             elseif sint == 12
#                 miu[5] = j
#             elseif sint == 6
#                 miu[6] = j
#             elseif sint == -4
#                 miu[7] = j
#             elseif sint == -2
#                 miu[8] = j
#             elseif sint == -10
#                 miu[9] = j
#             elseif sint == 8
#                 miu[10] = j
#             elseif sint == -12
#                 miu[11] = j
#             elseif sint == -6
#                 miu[12] = j
#             end
#         end
#     end
# # end
# #
# #         s = Integer.(sign.(rb) .* round.(rb/g))


#test for equilibrium: r0 =
# 1.476-1.476-0.000 2.952-1.476-1.476-1.476 1.476 0.000-2.952 1.476 1.476
# -0.852-2.557 1.704-0.000-0.852 2.557 0.852 2.557-1.704 0.000 0.852-2.557
# -2.410 0.000-2.410-0.000-2.410 0.000 2.410-0.000 2.410 0.000 2.410-0.000



# test_x = x0[100, :]
# test_r = r
#
# function test_potential(x)
#     V = 1.0/2.0*dot(x, d1_new * x)
#     return V
# end
#
#
#
# function test_force(x)
#     F = zeros(Float64, 3)
#     F = -1.0/2.0 * (d1_new + d1_new') * x
#     return F
# end
#
#
# a0 = test_force(test_r)
# q2 = grad(test_potential)
# a2 = -q2(r)
#
# q1 = x -> ForwardDiff.gradient(test_potential, x)
# a1 = -q1(r)
#
# a0./a2
