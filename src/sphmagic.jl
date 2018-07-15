function doSPHmagic(x, ρ, u, v, p, L, H)

    # 1 - DirectFind algorithm to find interacting particles
    (pairs, W, ∂x_W, W_0) = directfind(x)

    # 2 - SPH density approximation
    ρ = computedensity(pairs, W, W_0)

    # 3 - Equation of state for ideal gases
    p = (u .* ρ) * (GAMMA - 1)

    # 4,5 - Compute L operand
    (L, H) = computeoperands(ρ, p, v, pairs, ∂x_W)

    # 6 - Add artificial viscosity contribution
    (L, H) = correctwithartviscosity(x, ρ, p, v, pairs, ∂x_W, L, H)

    (x, ρ, u, v, p, L, H)
end

function directfind(x)
    pairs = []
    W = []
    ∂x_W = []

    for i in 1:PARTICLES_COUNT-1
        for j in i+1:PARTICLES_COUNT
            r = x[i] - x[j]
            if abs(r) < 2 * 2HSML
                push!(pairs, (i, j))

                (W_k, ∂x_W_k) = doublecosine(r, HSML)
                push!(W, W_k)
                push!(∂x_W, ∂x_W_k)
            end
        end
    end

    (W_0, _) = doublecosine(0, HSML)

    (pairs, W, ∂x_W, W_0)
end

function computedensity(pairs, W, W_0)
    # Considering the self effect
    ρ = ones(PARTICLES_COUNT) * PARTICLE_MASS * W_0

    # Computing for the neighbours
    for k in 1:length(pairs)
        (i, j) = pairs[k]
        ρ[i] += W[k] * PARTICLE_MASS
        ρ[j] += W[k] * PARTICLE_MASS
    end

    ρ
end

function computeoperands(ρ, p, v, pairs, ∂x_W)
    # Zero-initialization, since ∂x_W(0) = 0
    L = zeros(PARTICLES_COUNT)
    H = zeros(PARTICLES_COUNT)

    # Computing with the neighbours
    for k in 1:length(pairs)
        # Getting pairs
        (i, j) = pairs[k]

        # Calculating the constant term
        p_ij = p[i]/ρ[i]^2 + p[j]/ρ[j]^2

        # Wrapping for L and both particles, remembering ∂x_W_ij = -∂x_W_ji
        L[i] += ∂x_W[k] * p_ij * PARTICLE_MASS
        L[j] -= ∂x_W[k] * p_ij * PARTICLE_MASS

        # Wrapping for L and both particles, remembering ∂x_W_ij = -∂x_W_ji
        H[i] += ∂x_W[k] * 0.5 * (v[i] - v[j]) * p_ij * PARTICLE_MASS
        H[j] += ∂x_W[k] * 0.5 * (v[i] - v[j]) * p_ij * PARTICLE_MASS
    end

    (L, H)
end

function correctwithartviscosity(x, ρ, p, v, pairs, ∂x_W, L, H)
    # Constants
    qa = 1
    qb = 1
    etq = 0.1

    # Initialization of the corrections
    art_L = zeros(PARTICLES_COUNT)
    art_H = zeros(PARTICLES_COUNT)

    # Looping over pairs
    for k in 1:length(pairs)
        # Getting pairs
        (i, j) = pairs[k]

        x_ij = x[i] - x[j]
        v_ij = v[i] - v[j]

        xpsv = v_ij * x_ij
        xpsx = x_ij^2

        if (xpsv < 0)
            # Intermediary results
            muv = HSML * xpsv / (xpsx + HSML^2 * etq^2)
            mc = (sqrt(abs(GAMMA * p[i]/ρ[i])) + sqrt(abs(GAMMA * p[j]/ρ[j]))) / 2
            ρ_m = (ρ[i] + ρ[j]) / 2
            piv = (qb * muv - qb * mc) * muv/ρ_m
            hx = - piv * ∂x_W[k]

            # Adding to artificial acceleration contribution
            art_L[i] += PARTICLE_MASS * hx
            art_L[j] += PARTICLE_MASS * hx

            # Adding to artificial energy contribution
            art_H[i] += v_ij / 2 * PARTICLE_MASS * hx
            art_H[j] += v_ij / 2 * PARTICLE_MASS * hx
        end
    end

    # Applying the contribution
    L -= art_L
    H -= art_H

    (L, H)
end
