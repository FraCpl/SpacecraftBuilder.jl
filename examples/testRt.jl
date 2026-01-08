using LinearAlgebra
using JTools
using Quats
using SpacecraftBuilder
using GLMakie
using RayTracingEngine

function mainSc()
    cb = SpacecraftElement(;
        ID="CB",
        mass=6900.0,
        posEG_E=[25/1000, 3/1000, 1556/1000],
        geometry=Cuboid(2.5, 2.5, 1556*2/1000),
        inertiaG_E=[4285.0 -77 -311; -77 5070.0 102; -311 102 1795.0],
    )

    # Solar arrays
    θ = 90π/180
    mass = 212.0
    geom = Cuboid(6.648, 14.64, 0.015)
    inertiaG_E = [4730 -4.33 -0.026; -4.33 731 28; -0.026 28 5460]
    posEG_E = [0.0; 9.95; 0.02]

    sa1 = SpacecraftElement(;
        ID="SA (+Y)", mass=mass, geometry=geom, inertiaG_E=inertiaG_E, posEG_E=posEG_E, posOE_O=[0.0; 1.25; 1], R_OE=dcm_fromAxisAngle(2, θ)
    )

    sa2 = SpacecraftElement(;
        ID="SA (-Y)",
        mass=mass,
        geometry=geom,
        inertiaG_E=inertiaG_E,
        posEG_E=posEG_E,
        posOE_O=[0.0; -1.25; 1],
        R_OE=dcm_fromAxisAngle(3, π)*dcm_fromAxisAngle(2, -θ),
    )

    sc, model = build(cb, sa1, sa2; plotModel=true)

    # return rayTracingSrp(model, [-1.0; 0.0; 0.0]; Nrays=100_000)
end
mainSc()
