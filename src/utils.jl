@inline function rotateInertia(R_BA, J_A)
    J_B = similar(J_A)
    rotateInertia!(J_B, R_BA, J_A)
    return J_B
end

@inline function rotateInertia!(J_B, R_BA, J_A, tmp3x3=similar(J_A))
    mul!(tmp3x3, J_A, transpose(R_BA))
    mul!(J_B, R_BA, tmp3x3)
    return
end

@inline rotateModalMatrix(R_BA, L_A) = L_A*[R_BA' zeros(3, 3); zeros(3, 3) R_BA']         # L_B

# posGA = position of point A wrt element CoM (G).
# JG  = inertia at element CoM.
# posGA and JG must be written in the same reference frame (X).
#
# To get JG from JA, provide a negative mass!
# JG = translateInertia(JA,-mass,posGA)
@inline function translateInertia(JG_X, mass, posGA_X)
    JA_X = similar(JG_X)
    translateInertia!(JA_X, JG_X, mass, posGA_X)
    return JA_X
end

@inline function translateInertia!(JA_X, JG_X, mass, posGA_X)
    vvt!(JA_X, posGA_X) # faster than mul(JA_X, posGA_X, posGA_X')
    r2 = dot(posGA_X, posGA_X)
    JA_X[1, 1] -= r2
    JA_X[2, 2] -= r2
    JA_X[3, 3] -= r2
    @. JA_X = -mass*JA_X + JG_X
    # JA_X .*= -mass
    # JA_X .+= JG_X
    return
end

@inline function vvt!(M, v)
    @inbounds for i in 1:3, j in 1:3
        M[i, j] = v[i]*v[j]
    end
    return
end

@inline translateInertiaToCoM(JA_X, mass, posGA_X) = translateInertia(JA_X, -mass, posGA_X)     # JG_X

@inline function translateInertiaToCoM!(JG_X, JA_X, mass, posGA_X)
    translateInertia!(JG_X, JA_X, -mass, posGA_X)
    return
end

@inline translateModalMatrix(LA_X, posAB_X) = LA_X*[I crossmat(posAB_X); zeros(3, 3) I]     # LB_X

verifyInertia(J; verbose=true) = verifyInertia("_", J; verbose=verbose)

function verifyInertia(ID, J; verbose=true)
    if maximum(abs, J - J') > 1e-8
        if verbose
            @warn("Warning: the inertia matrix of $ID is not symmetric")
        end
        return false
    end

    eigsJ = eigvals(J)
    if any(.!isreal(eigsJ)) || any(eigsJ .< 0.0)
        if verbose
            @warn("Warning: the inertia matrix of $ID is not definite positive")
        end
        return false
    end

    if !(J[1, 1] ≤ J[2, 2] + J[3, 3] && J[2, 2] ≤ J[1, 1] + J[3, 3] && J[3, 3] ≤ J[1, 1] + J[2, 2])
        if verbose
            @warn("Warning: the diagonal terms of the inertia matrix of $ID are not physical")
        end
        return false
    end

    if !(abs(J[1, 2]) <= J[3, 3]/2 && abs(J[1, 3]) <= J[2, 2]/2 && abs(J[2, 3]) <= J[1, 1]/2)
        if verbose
            @warn("Warning: the non-diagonal terms of the inertia matrix of $ID are not physical")
        end
        return false
    end

    return true
end

# Caution: L and J shall be provided wrt the CoM
function verifyResidualMass(ID, m, JG, LG; verbose=true)

    # Necessary but not sufficient conditions
    r2 = [m; m; m; diag(JG)]
    l2 = sum(LG, dims=1)[:]
    del = r2 - l2 .< 0.0
    lbl = ["Ltx", "Lty", "Ltz", "Lθx", "Lθy", "Lθz"]
    for i in findall(del)
        if verbose
            @warn("Error with $(lbl[i])")
        end
    end
    if any(del); return false; end

    # Rigorous residual mass matrix check
    MR = [m*I zeros(3, 3); zeros(3, 3) JG] - LG'*LG
    if maximum(abs.(MR - MR')) > 1e-8
        if verbose
            @warn("Warning: the residual mass matrix of $ID is not symmetric")
        end
        return false
    end

    eigsJ = eigvals(MR)
    if any(.!isreal(eigsJ)) || any(eigsJ .< 0.0)
        if verbose
            @warn("Warning: the residual mass matrix of $ID is not definite positive ($(minimum(eigsJ)))")
        end
        return false
    end

    return true
end

randomInertia(Jnom, errPerc::Float64, rng) = randomInertia(Jnom, errPerc*randn(rng, 3))::Matrix{Float64}
randomInertia(Jnom, errPerc::Float64) = randomInertia(Jnom, errPerc*randn(3))::Matrix{Float64}
function randomInertia(Jnom::Matrix{Float64}, err::Vector{Float64})::Matrix{Float64}
    Jdiag, R = eigen(Jnom)
    A = (1.0 .+ err).*(sum(Jdiag)/2 .- Jdiag)
    return R*diagm(sum(A) .- A)*R'
end

function massMatrix2mci(MQ_B)
    if !isdiag(MQ_B[1:3, 1:3])
        @warn("Waring: mass not diagonal")
    end
    if norm(MQ_B - MQ_B') > 1e-9
        @warn("Mass matrix not symmetrical")
    end
    mass = MQ_B[1, 1]
    JQ_B = MQ_B[4:6, 4:6]
    verifyInertia(JQ_B)
    rQG_B = [-MQ_B[5, 3]; MQ_B[4, 3]; -MQ_B[4, 2]]./mass
    return mass, rQG_B, JQ_B
end
