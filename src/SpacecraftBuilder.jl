# Author: F. Capolupo
# European Space Agency, 2024
module SpacecraftBuilder

using JTools
using LinearAlgebra
using GLMakie
using GeometryBasics

export verifyInertia, randomInertia, translateInertia, rotateInertia
include("utils.jl")

export NoGeometry, PlanarPlate, Cuboid, Cylinder, Sphere
include("scGeometry.jl")

export SpacecraftElement
include("scElements.jl")

export build, getLTI, getMDK
include("scBuild.jl")

export rayTracingDrag, rayTracingSrp, rayTracingSurface
include("rayTracingUtils.jl")

end
