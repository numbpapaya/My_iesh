using My_iesh
using Test


@testset "My_iesh.jl" begin
    # Write your tests here.
    @test f(2, 1) == 7
    @test f(2, 3) == 13
    @test f(2, 3) == 1
end
