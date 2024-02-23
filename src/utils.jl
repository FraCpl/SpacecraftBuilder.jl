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

function verifyInertia(ID, J)
    if maximum(abs.(J - J')) > 1e-8
        println("Warning: the inertia matrix of $ID is not symmetric")
        return
    end

    eigsJ = eigvals(J)
    if any(.!isreal(eigsJ)) || any(eigsJ .< 0.0)
        println("Warning: the inertia matrix of $ID is not definite positive")
        return
    end

    if !(J[1, 1] ≤ J[2, 2] + J[3, 3] && J[2, 2] ≤ J[1, 1] + J[3, 3] && J[3, 3] ≤ J[1, 1] + J[2, 2])
        println("Warning: the diagonal terms of the inertia matrix of $ID are not physical")
        return
    end

    if !(abs(J[1, 2]) <= J[3, 3]/2 && abs(J[1, 3]) <= J[2, 2]/2 && abs(J[2, 3]) <= J[1, 1]/2)
        println("Warning: the non-diagonal terms of the inertia matrix of $ID are not physical")
        return
    end
end

function verifyResidualMass(ID, m, J, L)
    MR = [m*I zeros(3, 3); zeros(3, 3) J] - L'*L
    if maximum(abs.(MR - MR')) > 1e-8
        println("Warning: the residual mass matrix of $ID is not symmetric")
        return
    end

    eigsJ = eigvals(MR)
    if any(.!isreal(eigsJ)) || any(eigsJ .< 0.0)
        println("Warning: the residual mass matrix of $ID is not definite positive")
        return
    end
end
