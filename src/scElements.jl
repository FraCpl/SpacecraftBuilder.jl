
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
    freq::Vector{Float64}           # Input
    damp::Vector{Float64}           # Input
end

function SpacecraftElement(;
    ID::String="Element",
    mass::Float64=0.0,
    inertiaE_E::Matrix{Float64}=zeros(3, 3),        # [OPT1] Inertia
    inertiaG_E::Matrix{Float64}=zeros(3, 3),        # [OPT2] Inertia
    geometry::SpacecraftGeometry=NoGeometry(),
    R_OE::Matrix{Float64}=Matrix(1.0I, 3, 3),
    posEG_E::Vector{Float64}=zeros(3),
    posOE_O::Vector{Float64}=zeros(3),
    freq::Vector{Float64}=[NaN],
    damp::Vector{Float64}=[NaN],
    LE_E::Matrix{Float64}=zeros(length(freq), 6),      # [Translation, Rotation]
    LG_E::Matrix{Float64}=zeros(length(freq), 6),      # [Translation, Rotation]
)

    # CoM position
    posOG_O = posOE_O + R_OE*posEG_E

    # Compute inertias
    if any(inertiaE_E .!= 0.0)
        inertiaG_E .= translateInertiaToCoM(inertiaE_E, mass, posEG_E)
    elseif any(inertiaG_E .!= 0.0)
        inertiaE_E .= translateInertia(inertiaG_E, mass, posEG_E)
    else
        inertiaG_E .= elementInertiaG_E(geometry, mass)
        inertiaE_E .= translateInertia(inertiaG_E, mass, posEG_E)
    end
    verifyInertia(ID * " at E", inertiaE_E)
    verifyInertia(ID * " at G", inertiaG_E)
    inertiaG_O = rotateInertia(R_OE, inertiaG_E)
    inertiaO_O = translateInertia(inertiaG_O, mass, posOG_O)

    if isnan(freq[1])
        # Return rigid element
        return RigidElement(ID, geometry, mass, inertiaE_E, inertiaG_E, R_OE, posEG_E, posOE_O, posOG_O, inertiaO_O, inertiaG_O)
    end

    # Modal participation matrix
    if any(LE_E .!= 0.0)
        LG_E .= translateModalMatrix(LE_E, posEG_E)
    else
        LE_E .= translateModalMatrix(LG_E, -posEG_E)
    end
    LG_O = rotateModalMatrix(R_OE, LG_E)
    LO_O = translateModalMatrix(LG_O, -posOG_O)

    # Verify residual mass matrix
    verifyResidualMass(ID, mass, inertiaG_O, LG_O)

    # Sort things out
    id = sortperm(freq)

    # Return flexible element
    return FlexibleElement(
        ID,
        geometry,
        mass,
        inertiaE_E,
        inertiaG_E,
        R_OE,
        posEG_E,
        posOE_O,
        posOG_O,
        inertiaO_O,
        inertiaG_O,
        LO_O[id, :],
        LG_O[id, :],
        freq[id],
        damp[id],
    )
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
    #    println("freq:          $(round.(obj.freq, digits=2)) rad/s")
    #    println("damp:          $(round.(obj.damp*100, digits=2)) %")
    #    println("LO_O:       $(round.(obj.LO_O, digits=2))")
    #end
end
