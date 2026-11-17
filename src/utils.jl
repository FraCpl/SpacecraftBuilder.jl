@inline function invertInertiaCore(J)
    Jxx, Jxy, Jxz, ~, Jyy, Jyz, ~, ~, Jzz = J
    dj = 1/(Jxx*(Jyy*Jzz - (Jyz^2)) - Jxy*(Jxy*Jzz - Jxz*Jyz) + (Jxy*Jyz - Jxz*Jyy)*Jxz)

    # Compute inv(J)
    Ixx = (Jyy*Jzz - Jyz*Jyz)*dj
    Ixy = (-Jxy*Jzz + Jxz*Jyz)*dj
    Ixz = (Jxy*Jyz - Jxz*Jyy)*dj
    Iyy = (Jxx*Jzz - Jxz*Jxz)*dj
    Iyz = (-Jxx*Jyz + Jxy*Jxz)*dj
    Izz = (Jxx*Jyy - Jxy*Jxy)*dj
    return Ixx, Ixy, Ixz, Ixy, Iyy, Iyz, Ixz, Iyz, Izz
end

@inline function rotateInertiaCore(R_BA, J_A)
    # This implements J_B = R_BA*J_A*R_AB
    r11, r21, r31, r12, r22, r32, r13, r23, r33 = R_BA
    Jxx, Jxy, Jxz, ~, Jyy, Jyz, ~, ~, Jzz = J_A

    c1 = Jxx*r11 + Jxy*r12 + Jxz*r13
    c2 = Jxy*r11 + Jyy*r12 + Jyz*r13
    c3 = Jxz*r11 + Jyz*r12 + Jzz*r13
    c4 = Jxx*r21 + Jxy*r22 + Jxz*r23
    c5 = Jxy*r21 + Jyy*r22 + Jyz*r23
    c6 = Jxz*r21 + Jyz*r22 + Jzz*r23

    # Compute J_B
    Ixx = c1*r11 + c2*r12 + c3*r13
    Ixy = c1*r21 + c2*r22 + c3*r23
    Ixz = c1*r31 + c2*r32 + c3*r33
    Iyy = c4*r21 + c5*r22 + c6*r23
    Iyz = c4*r31 + c5*r32 + c6*r33
    Izz = (Jxx*r31 + Jxy*r32 + Jxz*r33)*r31 + (Jxy*r31 + Jyy*r32 + Jyz*r33)*r32 + (Jxz*r31 + Jyz*r32 + Jzz*r33)*r33
    return Ixx, Ixy, Ixz, Ixy, Iyy, Iyz, Ixz, Iyz, Izz
end

# This implements JA_X = JG_X - mass*crossMat(posGA_X)^2
@inline function translateInertiaCore(JG_X, mass, posGA_X)
    Jxx, Jxy, Jxz, ~, Jyy, Jyz, ~, ~, Jzz = JG_X
    x, y, z = posGA_X

    # Compute JA_X
    Ixx = Jxx + mass*(y*y + z*z)
    Ixy = Jxy - mass*x*y
    Ixz = Jxz - mass*x*z
    Iyy = Jyy + mass*(x*x + z*z)
    Iyz = Jyz - mass*y*z
    Izz = Jzz + mass*(x*x + y*y)
    return Ixx, Ixy, Ixz, Ixy, Iyy, Iyz, Ixz, Iyz, Izz
end

@inline invertInertia(J::Matrix3x3) = invertInertiaCore(J)
@inline invertInertia!(invJ, J) = toMatrix!(invJ, invertInertiaCore(J))

@inline rotateInertia(R_BA::Matrix3x3, J_A::Matrix3x3) = rotateInertiaCore(R_BA, J_A)
@inline rotateInertia(R_BA, J_A::AbstractMatrix) = toMatrix(rotateInertiaCore(R_BA, J_A))
@inline rotateInertia!(J_B, R_BA, J_A) = toMatrix!(J_B, rotateInertiaCore(R_BA, J_A))

@inline function q_rotateInertia!(J_B, q_BA, J_A)
    qs, qx, qy, qz = q_BA
    R_BA = q_toDcmCore(qs, qx, qy, qz)
    return rotateInertia!(J_B, R_BA, J_A)
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
@inline translateInertia(JG_X, mass, posGA_X) = toMatrix(translateInertiaCore(JG_X, mass, posGA_X))
@inline translateInertia!(JA_X, JG_X, mass, posGA_X) = toMatrix!(JA_X, translateInertiaCore(JG_X, mass, posGA_X))

@inline translateInertiaToCoM(JA_X, mass, posGA_X) = translateInertia(JA_X, -mass, posGA_X)     # JG_X
@inline translateInertiaToCoM!(JG_X, JA_X, mass, posGA_X) = translateInertia!(JG_X, JA_X, -mass, posGA_X)

@inline translateModalMatrix(LA_X, posAB_X) = LA_X*[I crossmat(posAB_X); zeros(3, 3) I]     # LB_X

verifyInertia(J; verbose=true) = verifyInertia("_", J; verbose=verbose)

function verifyInertia(ID, J; verbose=true)
    SYMMETRY_TOL = 1e-8
    Jxx, Jyx, Jzx, Jxy, Jyy, Jzy, Jxz, Jyz, Jzz = J

    val = abs(Jxy - Jyx)
    if val > SYMMETRY_TOL
        if verbose
            @warn("Warning: the inertia matrix of $ID is not symmetric: Jxy ≠ Jyx ($val)")
        end
        return false
    end

    val = abs(Jxz - Jzx)
    if val > SYMMETRY_TOL
        if verbose
            @warn("Warning: the inertia matrix of $ID is not symmetric: Jxz ≠ Jzx ($val)")
        end
        return false
    end

    val = abs(Jyz - Jzy)
    if val > SYMMETRY_TOL
        if verbose
            @warn("Warning: the inertia matrix of $ID is not symmetric: Jyz ≠ Jzy ($val)")
        end
        return false
    end

    for ej in eigvals(J)
        if ej < 0 || !isreal(ej)
            if verbose
                @warn("Warning: the inertia matrix of $ID is not definite positive ($ej)")
            end
            return false
        end
    end

    val1 = Jyy + Jzz - Jxx
    val2 = Jxx + Jzz - Jyy
    val3 = Jxx + Jyy - Jzz
    if !(val1 ≥ 0 && val2 ≥ 0 && val3 ≥ 0)
        if verbose
            @warn("Warning: the diagonal terms of the inertia matrix of $ID are not physical ($val1, $val2, $val3)")
        end
        return false
    end

    val1 = Jzz/2 - abs(Jxy)
    val2 = Jyy/2 - abs(Jxz)
    val3 = Jxx/2 - abs(Jyz)
    if !(val1 ≥ 0 && val2 ≥ 0 && val3 ≥ 0)
        if verbose
            @warn("Warning: the non-diagonal terms of the inertia matrix of $ID are not physical ($val1, $val2, $val3)")
        end
        return false
    end

    return true
end

# Caution: L and J shall be provided wrt the CoM
function verifyResidualMass(ID, m, JG, LG; verbose=true)

    # Necessary but not sufficient conditions
    r2 = [m; m; m; diag(JG)]
    l2 = sum(LG; dims=1)[:]
    del = r2 - l2 .< 0.0
    lbl = ["Ltx", "Lty", "Ltz", "Lθx", "Lθy", "Lθz"]
    for i in findall(del)
        if verbose
            @warn("Error with $(lbl[i])")
        end
    end
    if any(del)
        return false
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
            @warn("Warning: the residual mass matrix of $ID is not definite positive ($(minimum(eigsJ)))")
        end
        return false
    end

    return true
end

# Caution, when randomly sampling both the inertia and the CoM position, one must use the
# inertia at the CoG in the random sampling process, so that both JG and rOG can be sampled
# independently and still yield a physical configuration.
function randomInertia(JnomG, errPerc::Float64, rng)
    randomInertia(JnomG, errPerc*randn(rng, 3))::Matrix{Float64}
end

function randomInertia(JnomG, errPerc::Float64)
    randomInertia(JnomG, errPerc*randn(3))::Matrix{Float64}
end

function randomInertia(JnomG::Matrix{Float64}, err::Vector{Float64})::Matrix{Float64}
    Jdiag, R = eigen(JnomG)
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
