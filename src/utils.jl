@inline function invertInertia!(invJ, J)
    Jxx, Jxy, Jxz, ~, Jyy, Jyz, ~, ~, Jzz = J
    dj = 1/(Jxx*(Jyy*Jzz - (Jyz^2)) - Jxy*(Jxy*Jzz - Jxz*Jyz) + (Jxy*Jyz - Jxz*Jyy)*Jxz)
    invJ[1, 1] = (Jyy*Jzz - Jyz*Jyz)*dj
    invJ[1, 2] = (-Jxy*Jzz + Jxz*Jyz)*dj
    invJ[1, 3] = (Jxy*Jyz - Jxz*Jyy)*dj
    invJ[2, 2] = (Jxx*Jzz - Jxz*Jxz)*dj
    invJ[2, 3] = (-Jxx*Jyz + Jxy*Jxz)*dj
    invJ[3, 3] = (Jxx*Jyy - Jxy*Jxy)*dj

    invJ[2, 1] = invJ[1, 2]
    invJ[3, 1] = invJ[1, 3]
    invJ[3, 2] = invJ[2, 3]
    return
end

@inline function rotateInertia(R_BA, J_A)
    J_B = similar(J_A)
    rotateInertia!(J_B, R_BA, J_A)
    return J_B
end

@inline function rotateInertia!(J_B, R_BA, J_A)
    # This implements J_B = R_BA*J_A*R_AB
    r11, r21, r31, r12, r22, r32, r13, r23, r33 = R_BA
    Jxx, Jxy, Jxz, ~, Jyy, Jyz, ~, ~, Jzz = J_A

    c1 = Jxx*r11 + Jxy*r12 + Jxz*r13
    c2 = Jxy*r11 + Jyy*r12 + Jyz*r13
    c3 = Jxz*r11 + Jyz*r12 + Jzz*r13
    c4 = Jxx*r21 + Jxy*r22 + Jxz*r23
    c5 = Jxy*r21 + Jyy*r22 + Jyz*r23
    c6 = Jxz*r21 + Jyz*r22 + Jzz*r23

    J_B[1, 1] = c1*r11 + c2*r12 + c3*r13
    J_B[1, 2] = c1*r21 + c2*r22 + c3*r23
    J_B[1, 3] = c1*r31 + c2*r32 + c3*r33
    J_B[2, 2] = c4*r21 + c5*r22 + c6*r23
    J_B[2, 3] = c4*r31 + c5*r32 + c6*r33
    J_B[3, 3] =
        (Jxx*r31 + Jxy*r32 + Jxz*r33)*r31 +
        (Jxy*r31 + Jyy*r32 + Jyz*r33)*r32 +
        (Jxz*r31 + Jyz*r32 + Jzz*r33)*r33
    J_B[2, 1] = J_B[1, 2]
    J_B[3, 1] = J_B[1, 3]
    J_B[3, 2] = J_B[2, 3]
    return
end

# @inline function rotateInertia(R_BA::SMatrix{3, 3, T}, J_A::SMatrix{3, 3, T}) where T
#     return R_BA*J_A*transpose(R_BA)
# end

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
    # This implements JA_X = JG_X - mass*crossMat(posGA_X)^2
    Jxx, Jxy, Jxz, ~, Jyy, Jyz, ~, ~, Jzz = JG_X
    x, y, z = posGA_X

    JA_X[1, 1] = Jxx + mass*(y*y + z*z)
    JA_X[1, 2] = Jxy - mass*x*y
    JA_X[1, 3] = Jxz - mass*x*z
    JA_X[2, 2] = Jyy + mass*(x*x + z*z)
    JA_X[2, 3] = Jyz - mass*y*z
    JA_X[3, 3] = Jzz + mass*(x*x + y*y)
    JA_X[2, 1] = JA_X[1, 2]
    JA_X[3, 1] = JA_X[1, 3]
    JA_X[3, 2] = JA_X[2, 3]
    return
end

# @inline function translateInertia(JG_X::SMatrix{3, 3, T}, mass, posGA_X::SVector{3, T}) where T
#     return JG_X - mass*crossMatSq(posGA_X)
# end

@inline translateInertiaToCoM(JA_X, mass, posGA_X) = translateInertia(JA_X, -mass, posGA_X)     # JG_X

@inline function translateInertiaToCoM!(JG_X, JA_X, mass, posGA_X)
    translateInertia!(JG_X, JA_X, -mass, posGA_X)
    return
end

@inline translateModalMatrix(LA_X, posAB_X) = LA_X*[I crossmat(posAB_X); zeros(3, 3) I]     # LB_X

verifyInertia(J; verbose = true) = verifyInertia("_", J; verbose = verbose)

function verifyInertia(ID, J; verbose = true)
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

    if !(
        J[1, 1] ≤ J[2, 2] + J[3, 3] &&
        J[2, 2] ≤ J[1, 1] + J[3, 3] &&
        J[3, 3] ≤ J[1, 1] + J[2, 2]
    )
        if verbose
            @warn(
                "Warning: the diagonal terms of the inertia matrix of $ID are not physical"
            )
        end
        return false
    end

    if !(
        abs(J[1, 2]) <= J[3, 3]/2 && abs(J[1, 3]) <= J[2, 2]/2 && abs(J[2, 3]) <= J[1, 1]/2
    )
        if verbose
            @warn(
                "Warning: the non-diagonal terms of the inertia matrix of $ID are not physical"
            )
        end
        return false
    end

    return true
end

# Caution: L and J shall be provided wrt the CoM
function verifyResidualMass(ID, m, JG, LG; verbose = true)

    # Necessary but not sufficient conditions
    r2 = [m; m; m; diag(JG)]
    l2 = sum(LG, dims = 1)[:]
    del = r2 - l2 .< 0.0
    lbl = ["Ltx", "Lty", "Ltz", "Lθx", "Lθy", "Lθz"]
    for i in findall(del)
        if verbose
            @warn("Error with $(lbl[i])")
        end
    end
    if any(del)
        ;
        return false;
    end

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
            @warn(
                "Warning: the residual mass matrix of $ID is not definite positive ($(minimum(eigsJ)))"
            )
        end
        return false
    end

    return true
end

randomInertia(Jnom, errPerc::Float64, rng) =
    randomInertia(Jnom, errPerc*randn(rng, 3))::Matrix{Float64}
randomInertia(Jnom, errPerc::Float64) =
    randomInertia(Jnom, errPerc*randn(3))::Matrix{Float64}
function randomInertia(Jnom::Matrix{Float64}, err::Vector{Float64})::Matrix{Float64}
    Jdiag, R = eigen(Jnom)
    A = (1.0 .+ err) .* (sum(Jdiag)/2 .- Jdiag)
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
    rQG_B = [-MQ_B[5, 3]; MQ_B[4, 3]; -MQ_B[4, 2]] ./ mass
    return mass, rQG_B, JQ_B
end
