# Author: F. Capolupo
# European Space Agency, 2024
module SpacecraftBuilder

using JTools
using LinearAlgebra
using GLMakie
using GeometryBasics

export verifyInertia, randomInertia, translateInertia, rotateInertia, massMatrix2mci, translateInertiaToCoM
export translateInertia!, translateInertiaToCoM!, rotateInertia!
include("utils.jl")

export NoGeometry, PlanarPlate, Cuboid, Cylinder, Sphere
include("scGeometry.jl")

export SpacecraftElement
include("scElements.jl")

export build, getLTI, getMDK
include("scBuild.jl")

export rayTracingDrag, rayTracingSrp, rayTracingSurface, rayTracingHypersonicAero
include("rayTracingUtils.jl")

end
