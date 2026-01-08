function build(elements::Vararg{SpacecraftElement}; plotModel=false)
    if plotModel
        set_theme!(theme_fra(true))
        f = Figure(; size=(1500, 1000));
        display(f)
        p3d = LScene(f[1, 1]; show_axis=false)
        model = GeometryBasics.Mesh
        isFirst = true
    end

    # Init parameters
    ID = ""
    mass = 0.0
    posOG_O = zeros(3)
    inertiaO_O = zeros(3, 3)
    LO_O = zeros(1, 6)
    freq = [NaN]
    damp = [NaN]
    isFlex = false

    # Cycle through different elements
    for el in elements
        if el.mass > 0
            ID = ID*" + "*el.ID
            mass += el.mass
            posOG_O .+= el.mass*el.posOG_O
            inertiaO_O .+= el.inertiaO_O

            if el.isFlex
                LO_O = [LO_O; copy(el.LO_O)]
                freq = [freq; copy(el.freq)]
                damp = [damp; copy(el.damp)]

                if !isFlex
                    # Remove NaN values
                    LO_O = LO_O[2:end, :]
                    freq = freq[2:end]
                    damp = damp[2:end]
                    isFlex = true
                end
            end

            if plotModel
                if typeof(el.geometry) == Cuboid
                    if isFirst
                        model = cuboidModel(el.geometry.lx, el.geometry.ly, el.geometry.lz; pos_I=el.posOG_O, R_IB=el.R_OE)
                        isFirst = false
                    else
                        model = mergeMesh(model, cuboidModel(el.geometry.lx, el.geometry.ly, el.geometry.lz; pos_I=el.posOG_O, R_IB=el.R_OE))
                    end
                    #mesh!(p3d, cuboidModel(el.geometry.lx, el.geometry.ly, el.geometry.lz; pos_I=el.posOG_O, R_IB=el.R_OE); color=:lawngreen, alpha=0.5)
                end
                plotframe!(p3d, el.posOG_O, el.R_OE, 2.0)
                scatter!(p3d, el.posOE_O[1], el.posOE_O[2], el.posOE_O[3]; markersize=10)   # This gets hidden by the mesh!
            end
        end
    end

    # Normalize CoM position wrt total mass
    posOG_O ./= mass

    if plotModel
        mesh!(p3d, model; color=:lawngreen, alpha=0.5)
        plotframe!(p3d, posOG_O, Matrix(1.0I, 3, 3), 3.0)
    end

    # Return assembled spacecraft element
    sc = SpacecraftElement(; ID=ID[4:end], mass=mass, inertiaE_E=inertiaO_O, posEG_E=posOG_O, freq=freq, damp=damp, LE_E=LO_O)

    if plotModel
        return sc, model
    end
    return sc
end

function getLTI(elements::Vector; attitudeOnly=false)
    return getLTI(build(elements); attitudeOnly=attitudeOnly)
end

function getLTI(sc::SpacecraftElement; attitudeOnly=false)
    M, D, K = getMDK(sc)
    nf = sc.isFlex ? length(sc.freq) : 0
    if attitudeOnly
        nu = 3
        Bu = [zeros(3, nu); I; zeros(nf, nu)]
        Cy = [zeros(3, 3) I zeros(3, nf)]
    else
        nu = 6
        Bu = [I; zeros(nf, nu)]
        Cy = [I zeros(6, nf)]
    end
    A = [zeros(6 + nf, 6 + nf) I; -M\K -M\D]
    B = [zeros(6 + nf, nu); M\Bu]
    C = [Cy zeros(size(Cy))]
    return A, B, C, 0
end

function getMDK(sc::SpacecraftElement)
    if !sc.isFlex
        return [sc.mass*I zeros(3, 3); zeros(3, 3) sc.inertiaG_O], zeros(6, 6), zeros(6, 6)
    end
    M = [[sc.mass*I zeros(3, 3); zeros(3, 3) sc.inertiaG_O] sc.LG_O'; sc.LG_O I]
    D = diagm([zeros(6); 2sc.damp .* sc.freq])
    K = diagm([zeros(6); sc.freq .^ 2])
    return M, D, K
end
