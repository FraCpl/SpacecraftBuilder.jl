using SpacecraftBuilder
using Test
using Quaternions
using LinearAlgebra

@testset "SpacecraftBuilder.jl" begin
    # Write your tests here.
    JG = randn(3, 3)
    m = rand()
    r = randn(3)
    R = randn(3, 3)
    @test norm(translateInertia(JG, m, r) - (JG - m*crossMat(r)^2)) < 1e-12
    @test norm(rotateInertia(R, JG) - R*JG*R') < 1e-12
    @test norm(translateInertiaToCoM(JG, m, r) - (JG + m*crossMat(r)^2)) < 1e-12
end
