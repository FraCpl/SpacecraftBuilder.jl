module SpacecraftBuilder

using JTools
using LinearAlgebra
using GLMakie
using ControlSystems

export verifyInertia
include("utils.jl")

export NoGeometry, PlanarPlate, Cuboid, Cylinder, Sphere
include("scGeometry.jl")

export RigidElement, FlexibleElement
include("scElements.jl")

export build, buildss
include("scBuild.jl")

end
