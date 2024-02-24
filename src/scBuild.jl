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

    # Init parameters
    ID = ""
    mass = 0.0
    posOG_O = zeros(3)
    inertiaO_O = zeros(3, 3)
    LO_O = zeros(1, 6)
    ω = [NaN]
    ξ = [NaN]
    hasFlex = false

    # Cycle through different elements
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

    # Normalize CoM position wrt total mass
    posOG_O ./= mass

    if plotModel
        plotframe!(p3d, posOG_O, Matrix(1.0I, 3, 3), 3.0)
        display(f)
    end

    # Return assembled spacecraft element
    return SpacecraftElement(ID=ID[4:end], mass=mass, inertiaE_E=inertiaO_O, posEG_E=posOG_O,
        ω=ω, ξ=ξ, LE_E=LO_O)
end

function buildss(elements::Vector; attitudeOnly=false)
    return buildss(build(elements); attitudeOnly=attitudeOnly)
end

function buildss(sc::SpacecraftElement; attitudeOnly=false)
    nω = length(sc.ω)
    M, D, K = getMDK(sc)
    if attitudeOnly
        nu = 3
        Bu = [zeros(3, nu); I; zeros(nω, nu)]
        Cy = [zeros(3, 3) I zeros(3, nω)]
    else
        nu = 6
        Bu = [I; zeros(nω, nu)]
        Cy = [I zeros(6, nω)]
    end
    A = [zeros(6 + nω, 6 + nω) I; -M\K -M\D]
    B = [zeros(6 + nω, nu); M\Bu]
    C = [Cy zeros(size(Cy))]
    return ss(A, B, C, 0)
end

function getMDK(sc::SpacecraftElement)
    if typeof(sc) == RigidElement
        return [sc.mass*I zeros(3, 3); zeros(3, 3) sc.inertiaG_O], zeros(6, 6), zeros(6, 6)
    end
    M = [[sc.mass*I zeros(3, 3); zeros(3, 3) sc.inertiaG_O] sc.LG_O'; sc.LG_O I]
    D = diagm([zeros(6); 2sc.ξ.*sc.ω])
    K = diagm([zeros(6); sc.ω.^2])
    return M, D, K
end
