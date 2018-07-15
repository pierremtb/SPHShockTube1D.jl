include("sphmagic.jl")

function computeSimulationStep(step, x, ρ, u, v, p, L, H)
    # First Leap-Frog integration
    if step > 1
        # Keep track of previous energy and velocity
        previous_u = copy(u)
        previous_v = copy(v)

        # Compute half-step integration
        u += TIME_DELTA/2 * H
        v += TIME_DELTA/2 * L
    end

    (x, ρ, u, v, p, L, H) = doSPHmagic(x, ρ, u, v, p, L, H)

    # Second Leap-Frog integration
    if step > 1
        # Compute half-step integration with previous value
        u = previous_u + TIME_DELTA * H
        v = previous_v + TIME_DELTA * L
    else
        # Compute first half-step integration
        u += TIME_DELTA * H
        v += TIME_DELTA * L
    end
    # Compute regular integration
    x += TIME_DELTA * v

    (x, ρ, u, v, p, L, H)
end
