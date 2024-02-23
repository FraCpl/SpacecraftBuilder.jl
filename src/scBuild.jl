function initElement(;
        ID::String="Element",
        mass::Float64=0.0,
        inertiaE_E::Matrix{Float64}=zeros(3, 3),        # [OPT1] Inertia
        inertiaG_E::Matrix{Float64}=zeros(3, 3),        # [OPT2] Inertia
        geometry::SpacecraftGeometry=NoGeometry(),
        R_OE::Matrix{Float64}=Matrix(1.0I, 3, 3),
        posEG_E::Vector{Float64}=zeros(3),
        posOE_O::Vector{Float64}=zeros(3),
        ω::Vector{Float64}=[NaN],
        ξ::Vector{Float64}=[NaN],
        LE_E::Matrix{Float64}=zeros(length(ω), 6),      # [OPT1] Modal participation matrix = [Translation, Rotation]
        LG_E::Matrix{Float64}=zeros(length(ω), 6),      # [OPT2] Modal participation matrix = [Translation, Rotation]
    )

    # CoM position
    posOG_O = posOE_O + R_OE*posEG_E

    # Compute inertias
    if any(inertiaE_E .!= 0.0)
        verifyInertia(ID, inertiaE_E)
        inertiaG_E .= translateInertia(inertiaE_E, -mass, posEG_E)
    elseif any(inertiaG_E .!= 0.0)
        verifyInertia(ID, inertiaG_E)
        inertiaE_E .= translateInertia(inertiaG_E, +mass, posEG_E)
    else
        inertiaG_E .= elementInertiaG_E(geometry, mass)
        inertiaE_E .= translateInertia(inertiaG_E, +mass, posEG_E)
    end
    inertiaG_O = rotateInertia(R_OE, inertiaG_E)
    inertiaO_O = translateInertia(inertiaG_O, mass, posOG_O)

    if isnan(ω[1])
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
    id = sortperm(ω)

    # Return flexible element
    return FlexibleElement(ID, geometry, mass, inertiaE_E, inertiaG_E, R_OE, posEG_E, posOE_O, posOG_O, inertiaO_O, inertiaG_O, LO_O[id, :], LG_O[id, :], ω[id], ξ[id])
end


#=
# from O to Q (=new O)
function transformElement!(el::SpacecraftElement, posQO_Q=zeros(3), R_QO=I)
el.R_OE .= R_QO*el.R_OE
el.posOE_O .= posQO_Q + R_QO*el.posOE_O
el.posOG_O = posQO_Q + R_QO*el.posOG_O
el.inertiaG_O = rotateInertia(R_QO, el.inertiaG_O)
el.inertiaO_O = translateInertia(el.inertiaG_O, el.mass, el.posOG_O)
if typeof(el) == FlexibleElement
    el.LG_O = rotateModalMatrix(R_QO, el.LG_O)
    el.LO_O = translateModalMatrix(el.LG_O, -el.posOG_O)
end
end

function transformElement(el::SpacecraftElement, posQO_Q=zeros(3), R_QO=I)
el2 = deepcopy(el)
transformElement!(el2, posQO_Q, R_QO)
return el2
end
=#

function build(elements::Vector; plotModel=false)

    if plotModel
        set_theme!(theme_fra(true))
        f = Figure(size=(1500, 1000))
        p3d = LScene(f[1, 1]; show_axis=false)
    end

    ID = ""
    mass = 0.0
    posOG_O = zeros(3)
    inertiaO_O = zeros(3, 3)
    LO_O = zeros(1, 6)
    ω = [NaN]
    ξ = [NaN]
    hasFlex = false
    for el in elements
        ID = ID*" + "*el.ID
        mass += el.mass
        posOG_O .+= el.mass*el.posOG_O
        inertiaO_O .+= el.inertiaO_O

        if typeof(el) == FlexibleElement
            LO_O = [LO_O; copy(el.LO_O)]
            ω = [ω; copy(el.ω)]
            ξ = [ξ; copy(el.ξ)]
            if !hasFlex
                # Remove NaN values
                LO_O = LO_O[2:end, :]
                ω = ω[2:end]
                ξ = ξ[2:end]
                hasFlex = true
            end
        end

        if plotModel
            if typeof(el.geometry) == Cuboid
                mesh!(p3d, cuboidModel(el.geometry.lx, el.geometry.ly, el.geometry.lz; pos_I=el.posOG_O, R_IB=el.R_OE); color=:lawngreen, alpha=0.5)
            end
            plotframe!(p3d, el.posOG_O, el.R_OE, 2.0)
            scatter!(p3d, el.posOE_O[1], el.posOE_O[2], el.posOE_O[3]; markersize=10)   # This gets hidden by the mesh!
        end
    end
    posOG_O ./= mass

    if plotModel
        plotframe!(p3d, posOG_O, Matrix(1.0I, 3, 3), 3.0)
        display(f)
    end

    return initElement(ID=ID[4:end], mass=mass, inertiaE_E=inertiaO_O, posEG_E=posOG_O, ω=ω, ξ=ξ, LE_E=LO_O)
end


function buildss(sc::SpacecraftElement; attitudeOnly=false)
    nω = length(sc.ω)
    if attitudeOnly
        nx = nω + 3
        nu = 3
        M = [sc.inertiaG_O sc.LG_O[:, 4:6]'; sc.LG_O[:, 4:6] I]
    else
        nx = nω + 6
        nu = 6
        M = [[sc.mass*I zeros(3, 3); zeros(3, 3) sc.inertiaG_O] sc.LG_O'; sc.LG_O I]
    end
    D = diagm([zeros(nx - nω); 2sc.ξ.*sc.ω])
    K = diagm([zeros(nx - nω); sc.ω.^2])
    Bu = [I; zeros(nω, nu)]
    Cy = [I zeros(nx - nω, nω)]
    A = [zeros(nx, nx) I; -M\K -M\D]
    B = [zeros(nx, nu); M\Bu]
    C = [Cy zeros(size(Cy))]
    return ss(A, B, C, 0)
end
