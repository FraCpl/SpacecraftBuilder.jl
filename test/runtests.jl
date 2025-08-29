using SpacecraftBuilder
using Test
using Quaternions
using LinearAlgebra

@testset "SpacecraftBuilder.jl" begin
    JG = randn(3, 3)
    JG = JG + JG'
    m = rand()
    r = randn(3)
    R = randn(3, 3)
    invJ = zeros(3, 3)
    invertInertia!(invJ, JG)
    @test norm(translateInertia(JG, m, r) - (JG - m*crossMat(r)^2)) < 1e-12
    @test norm(rotateInertia(R, JG) - R*JG*R') < 1e-12
    @test norm(translateInertiaToCoM(JG, m, r) - (JG + m*crossMat(r)^2)) < 1e-12
    @test norm(invJ - inv(JG)) < 1e-12
end
