include("sphmagic.jl")

function computeSimulationStep(step, x, ρ, u, v, p, L, H)
    # First Leap-Frog integration
    if step > 1
        # Keep track of previous energy and velocity
        previous_u = copy(u)
        previous_v = copy(v)

        # Compute half-step integration
        for i in FIXED_SIDES_PARTICLES_COUNT:PARTICLES_COUNT-FIXED_SIDES_PARTICLES_COUNT
            u[i] += TIME_DELTA/2 * H[i]
            v[i] += TIME_DELTA/2 * L[i]
        end
    end

    (x, ρ, u, v, p, L, H) = doSPHmagic(x, ρ, u, v, p, L, H)

    # Second Leap-Frog integration
    for i in FIXED_SIDES_PARTICLES_COUNT:PARTICLES_COUNT-FIXED_SIDES_PARTICLES_COUNT
        if step > 1
            # Compute half-step integration with previous value
            u[i] = previous_u[i] + TIME_DELTA * H[i]
            v[i] = previous_v[i] + TIME_DELTA * L[i]
        else
            # Compute first half-step integration
            u[i] += TIME_DELTA * H[i]
            v[i] += TIME_DELTA * L[i]
        end
        # Compute regular integration
        x[i] += TIME_DELTA * v[i]
    end

    (x, ρ, u, v, p, L, H)
end
