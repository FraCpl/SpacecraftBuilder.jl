
abstract type SpacecraftGeometry end

mutable struct NoGeometry <: SpacecraftGeometry
end

mutable struct Cuboid <: SpacecraftGeometry
    lx::Float64
    ly::Float64
    lz::Float64
end

mutable struct PlanarPlate <: SpacecraftGeometry
    lx::Float64
    ly::Float64
end

mutable struct Cylinder <: SpacecraftGeometry
    r::Float64
    h::Float64
end

mutable struct Sphere <: SpacecraftGeometry
    r::Float64
end

elementInertiaG_E(el::Cuboid, mass) = mass/12*[el.ly^2 + el.lz^2 0.0 0.0; 0.0 el.lx^2 + el.lz^2 0.0; 0.0 0.0 el.ly^2 + el.lx^2]
elementInertiaG_E(el::PlanarPlate, mass) = mass/12*[el.ly^2 0.0 0.0; 0.0 el.lx^2 0.0; 0.0 0.0 el.ly^2 + el.lx^2]
elementInertiaG_E(el::Cylinder, mass) = mass*[3/12*el.r^2 + 1/12*el.h^2 0.0 0.0; 0.0 3/12*el.r^2 + 1/12*el.h^2 0.0; 0.0 0.0 1/2*el.r^2]
elementInertiaG_E(el::Sphere, mass) = mass*2/5*el.r^2*Matrix(1.0I, 3, 3)
