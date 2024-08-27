mutable struct Ray
    origin::Vector
    dir::Vector
    t::Real
    idxFace::Int
end

Ray(origin, dir) = Ray(origin, normalize(dir), Inf, 0)

const EPS_PARALLEL = 1e-8
const ZERO_WITH_TOL = -1e-8
const ONE_WITH_TOL = 1.0 + 1e-8
const EPS_MIN_DIST = 1e-4
const RAY_SHIFT = 1e-2
const DEFAULT_NRAYS = 1_000_000

@views function intersect!(ray::Ray, tri, idx::Int)
    E12 = tri[2] - tri[1]
    E13 = tri[3] - tri[1]
    P = ray.dir × E13
    detM = dot(P, E12)
    if abs(detM) ≤ EPS_PARALLEL
        return
    end

    T = ray.origin - tri[1]
    u = dot(P, T)/detM
    if u < ZERO_WITH_TOL || u > ONE_WITH_TOL
        return
    end

    Q = T × E12
    v = dot(Q, ray.dir)/detM
    if v < ZERO_WITH_TOL || u + v > ONE_WITH_TOL
        return
    end

    t = dot(E13, Q)/detM
    if t > EPS_MIN_DIST && t < ray.t
        ray.t = t
        ray.idxFace = idx
    end
end

@views function intersect!(ray::Ray, bboxMin, bboxMax)
    t1 = (bboxMin - ray.origin)./ray.dir
    t2 = (bboxMax - ray.origin)./ray.dir
    tmin = min(t1[1], t2[1])
    tmax = max(t1[1], t2[1])

    tmin = max(tmin, min(t1[2], t2[2]))
    tmax = min(tmax, max(t1[2], t2[2]))

    tmin = max(0.0, max(tmin, min(t1[3], t2[3])))
    tmax = min(tmax, max(t1[3], t2[3]))
    return tmax ≥ tmin
end

# Qdyn = 1/2ρV^2
function rayTracingDrag(model::GeometryBasics.Mesh, dirVel::Vector{Float64}, Qdyn=1.0, Cx=1.0; Nrays::Int=DEFAULT_NRAYS)
    S, posCoP = rayTracingSrp(model, dirVel; Nrays=Nrays, Nrec=1, mode=:drag)
    forceDrag = -Qdyn*Cx*S*normalize(dirVel)
    torqueDrag = posCoP × forceDrag
    return forceDrag, torqueDrag
end

# Newton theory for hypersonic flow
# dirVel is the direction of the velocity of the body (and not the direction of the flow as`
# seen from the body).
#
# The output torque is computed with respect to the origin of the 3D model provided as input,
# or with respect to the point specified by the optional argument 'posRef'.
#
# To compute drag, lift, and torque coefficients instead of drag, lift and torque, just set
# the dynamic pressure to 1.0, and divide the outputs by the reference surface of the body
# (to be chosen by the user).
#
# References:
# [1] Anderson, Fundamentals of Aerodynamics
# [2] Heybey, Newtonian Aerodynamics for General Body Shapes with Several Applications,
#     https://ntrs.nasa.gov/api/citations/19660012440/downloads/19660012440.pdf
function rayTracingHypersonicAero(model::GeometryBasics.Mesh, dirVel::Vector{Float64}, Qdyn=1.0, CpMax=1.84; Nrays::Int=DEFAULT_NRAYS, posRef=zeros(3))
    dvel = normalize(dirVel)
    force, torque = rayTracingSrp(model, dvel; Nrays=Nrays, Nrec=1, mode=:hyper)
    force *= CpMax*Qdyn
    torque *= CpMax*Qdyn
    torque -= posRef × force
    drag = dot(force, -dvel).*(-dvel)
    lift = force - drag
    return drag, lift, torque
end

function rayTracingSurface(model::GeometryBasics.Mesh, dirObs::Vector{Float64}; Nrays::Int=DEFAULT_NRAYS)
    return rayTracingSrp(model, dirObs; Nrays=Nrays, Nrec=1, mode=:surf)[1]
end

@views function rayTracingSrp(model::GeometryBasics.Mesh, dirSun::Vector{Float64}, Psrp=1.0; Nrays::Int=DEFAULT_NRAYS, Nrec=3, mode=:srp)

    anyIntersection = !(mode == :srp || mode == :hyper)
    α = 0.7; rd = 0.1; rs = 0.2

    # Model box
    bboxMin = [minimum(getindex.(model.position, i)) for i in 1:3]
    bboxMax = [maximum(getindex.(model.position, i)) for i in 1:3]
    X0 = (bboxMin + bboxMax)/2
    R = maximum(bboxMax - X0)

    # Creates rays matrix
    coords = 1.1*R*range(-1, 1, round(Int, √Nrays))
    dirRay = -normalize(dirSun)
    R_IS = genOrthogonalAxes(dirSun)
    Ap = (coords[2] - coords[1])^2

    rays = Ray[]
    for x in coords, y in coords
        push!(rays, Ray(X0 + R_IS*[x; y; 3R], dirRay))
    end

    # Intersect model
    surf = 0.0; posCoP = zeros(3)
    forceSrp = zeros(3); torqueSrp = zeros(3)
    @inbounds for ray in rays
        # Check if the ray intersects the object's bounding box
        if intersect!(ray, bboxMin, bboxMax)
            psrp = Psrp
            # Initialize ray recursion
            @inbounds for n in 1:Nrec
                # Intersect ray with entire model (closest intersection)
                @inbounds for k in eachindex(model)
                    intersect!(ray, model[k], k)
                    if anyIntersection && ray.t < Inf; break; end
                end

                # Check if any intersection has been found
                if ray.t == Inf
                    break
                elseif mode == :srp
                    # Compute SRP
                    N, cosθ, newDir = rayReflection(model[ray.idxFace], ray)
                    frc, trq, psrp, posHit = srpFormula(ray, psrp, Ap, N, α, rd, rs, cosθ)
                    forceSrp += frc
                    torqueSrp += trq

                    # Update ray for next reflection
                    if n < Nrec
                        ray.t = Inf
                        ray.dir = newDir
                        ray.origin = posHit + RAY_SHIFT*newDir
                    end
                elseif mode == :hyper
                    # Newton theory for hypersonic aerodynamics
                    N, cosθ, ~ = rayReflection(model[ray.idxFace], ray)     # N is the outward normal
                    frc = -N*cosθ*Ap                                        # sinα = sin(π/2 - θ) = cosθ, Ap = dS*cosθ = dS*sinα,
                    posHit = ray.origin + ray.t*ray.dir
                    torqueSrp += posHit × frc
                    forceSrp += frc
                else
                    # Compute surface and center of pressure
                    surf += Ap
                    posCoP .+= ray.origin*Ap
                end
            end
        end
    end

    # Output
    if mode == :srp
        return forceSrp, torqueSrp
    elseif mode == :hyper
        return forceSrp, torqueSrp
    else
        return surf, posCoP./surf
    end
end

function rayReflection(face, ray)
    E1 = face[2] - face[1]
    E2 = face[3] - face[2]
    N = normalize(E1 × E2)
    cosθ = -N'*ray.dir

    # Fix normal direction of current face on the basis of ray
    # intersection result (the dot product between minus ray direction
    # and face normal must be positive)
    if cosθ < 0.0
        N = -N
        cosθ = -cosθ
    end

    return N, cosθ, ray.dir + 2cosθ*N
end

function srpFormula(ray, Psrp, Ap, N, α, rd, rs, cosθ)

    # Compute forces generated by each single ray intersecting the object
    forceSrp = -Psrp*Ap*(-(α + rd)*ray.dir + 2(rs*cosθ + rd/3)*N)

    # Compute total torque wrt the origin of the model
    # (i.e., origin of the coordinates of the vertex of the model)
    posHit = ray.origin + ray.t*ray.dir
    torqueSrp = posHit × forceSrp

    # Update SRP flux for next batch of reflected rays
    return forceSrp, torqueSrp, Psrp*rs, posHit
end

function genOrthogonalAxes(zB_A)
    nu = norm(zB_A)
    if nu > 0.0
        zB_A = zB_A./nu
        VB_A = [zB_A × [i==1; i==2; i==3] for i in 1:3]
        vNorm = norm.(VB_A)
        iMax = findfirst(vNorm .≥ 0.2)
        vB_A = VB_A[iMax]./vNorm[iMax]
        return [(vB_A × zB_A) vB_A zB_A]
    end
    return Matrix(1.0I, 3, 3)
end
