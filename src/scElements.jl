
abstract type SpacecraftElement end

mutable struct RigidElement <: SpacecraftElement
    ID::String
    geometry::SpacecraftGeometry
    mass::Float64                # Input
    inertiaE_E::Matrix{Float64}  # Possible input
    inertiaG_E::Matrix{Float64}  # Possible input
    R_OE::Matrix{Float64}        # Input
    posEG_E::Vector{Float64}     # Input
    posOE_O::Vector{Float64}     # Input
    posOG_O::Vector{Float64}     # Automatically computed
    inertiaO_O::Matrix{Float64}  # Automatically computed
    inertiaG_O::Matrix{Float64}  # Automatically computed
end

mutable struct FlexibleElement <: SpacecraftElement
    ID::String
    geometry::SpacecraftGeometry
    mass::Float64                   # Input
    inertiaE_E::Matrix{Float64}     # Possible input
    inertiaG_E::Matrix{Float64}     # Possible input
    R_OE::Matrix{Float64}           # Input
    posEG_E::Vector{Float64}        # Input
    posOE_O::Vector{Float64}        # Input
    posOG_O::Vector{Float64}        # Automatically computed
    inertiaO_O::Matrix{Float64}     # Automatically computed
    inertiaG_O::Matrix{Float64}     # Automatically computed
    LO_O::Matrix{Float64}           # Automatically computed
    LG_O::Matrix{Float64}           # Automatically computed
    ω::Vector{Float64}              # Input
    ξ::Vector{Float64}              # Input
end

function RigidElement(;
        ID::String="R-Element",
        mass::Float64=0.0,
        inertiaE_E::Matrix{Float64}=zeros(3, 3),        # [OPT1] Inertia
        inertiaG_E::Matrix{Float64}=zeros(3, 3),        # [OPT2] Inertia
        geometry::SpacecraftGeometry=NoGeometry(),
        R_OE::Matrix{Float64}=Matrix(1.0I, 3, 3),
        posEG_E::Vector{Float64}=zeros(3),
        posOE_O::Vector{Float64}=zeros(3)
    )

    return initElement(ID=ID, geometry=geometry, mass=mass, inertiaE_E=inertiaE_E,
        inertiaG_E=inertiaG_E, R_OE=R_OE, posEG_E=posEG_E, posOE_O=posOE_O)
end

function FlexibleElement(;
        ID::String="F-Element",
        mass::Float64=0.0,
        inertiaE_E::Matrix{Float64}=zeros(3, 3),        # [OPT1] Inertia
        inertiaG_E::Matrix{Float64}=zeros(3, 3),        # [OPT2] Inertia
        geometry::SpacecraftGeometry=NoGeometry(),
        R_OE::Matrix{Float64}=Matrix(1.0I, 3, 3),
        posEG_E::Vector{Float64}=zeros(3),
        posOE_O::Vector{Float64}=zeros(3),
        ω::Vector{Float64}=[NaN],
        ξ::Vector{Float64}=[NaN],
        LE_E::Matrix{Float64}=zeros(length(ω), 6),      # [Translation, Rotation]
        LG_E::Matrix{Float64}=zeros(length(ω), 6),      # [Translation, Rotation]
    )

    return initElement(ID=ID, geometry=geometry, mass=mass, inertiaE_E=inertiaE_E,
        inertiaG_E=inertiaG_E, R_OE=R_OE, posEG_E=posEG_E, posOE_O=posOE_O, ω=ω, ξ=ξ,
        LE_E=LE_E, LG_E=LG_E)
end

function Base.show(io::IO, obj::SpacecraftElement)
    println("ID:         $(obj.ID)")
    println("type:       $(typeof(obj))")
    println("mass:       $(round(obj.mass, digits=2)) kg")
    println("posOG_O:    $(round.(obj.posOG_O, digits=3)) m")
    println("inertia at CoM:")
    println("  Ixx, Iyy, Izz: $(round.(diag(obj.inertiaG_O), digits=1)) kg m²")
    println("  Ixy, Ixz, Iyz: $(round.([obj.inertiaG_O[1, 2], obj.inertiaG_O[1, 3], obj.inertiaG_O[2, 3]], digits=1)) kg m²")
    #if typeof(obj) == FlexibleElement
    #    println("ω:          $(round.(obj.ω, digits=2)) rad/s")
    #    println("ξ:          $(round.(obj.ξ*100, digits=2)) %")
    #    println("LO_O:       $(round.(obj.LO_O, digits=2))")
    #end
end
