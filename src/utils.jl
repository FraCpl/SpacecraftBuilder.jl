rotateInertia(R_BA, J_A) = R_BA*J_A*R_BA'  # J_B
rotateModalMatrix(R_BA, L_A) = L_A*[R_BA' zeros(3, 3); zeros(3, 3) R_BA']         # L_B

# posGA = position of point A wrt element CoM (G).
# JG  = inertia at element CoM.
# posGA and JG must be written in the same reference frame (X).
#
# To get JG from JA, provide a negative mass!
# JG = translateInertia(JA,-mass,posGA)
translateInertia(JG_X, mass, posGA_X) = JG_X - mass*crossmat(posGA_X)^2    # JA_X
translateModalMatrix(LA_X, posAB_X) = LA_X*[I crossmat(posAB_X); zeros(3, 3) I]     # LB_X

verifyInertia(J; verbose=true) = verifyInertia("_", J; verbose=verbose)

function verifyInertia(ID, J; verbose=true)
    if maximum(abs.(J - J')) > 1e-8
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

function verifyResidualMass(ID, m, J, L; verbose=true)
    MR = [m*I zeros(3, 3); zeros(3, 3) J] - L'*L
    if maximum(abs.(MR - MR')) > 1e-8
        if verbose
            @warn("Warning: the residual mass matrix of $ID is not symmetric")
        end
        return false
    end

    eigsJ = eigvals(MR)
    if any(.!isreal(eigsJ)) || any(eigsJ .< 0.0)
        if verbose
            @warn("Warning: the residual mass matrix of $ID is not definite positive")
        end
        return false
    end

    return true
end

randomInertia(Jnom, errPerc::Float64) = randomInertia(Jnom, errPerc*randn(3))::Matrix{Float64}
function randomInertia(Jnom::Matrix{Float64}, err::Vector{Float64})::Matrix{Float64}
    Jdiag, R = eigen(Jnom)
    A = (1.0 .+ err).*(sum(Jdiag)/2 .- Jdiag)
    return R*diagm(sum(A) .- A)*R'
end
