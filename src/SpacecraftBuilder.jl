module SpacecraftBuilder

using JTools
using LinearAlgebra
using GLMakie

export verifyInertia
include("utils.jl")

export NoGeometry, PlanarPlate, Cuboid, Cylinder, Sphere
include("scGeometry.jl")

export SpacecraftElement
include("scElements.jl")

export build, getLTI, getMDK
include("scBuild.jl")

end
