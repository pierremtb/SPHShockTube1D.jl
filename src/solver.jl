using Base
using Plots

include("initialconditions.jl")
include("kernels.jl")
include("output.jl")
include("simulationstep.jl")

function runsimulation()
    # Initial quantities
    (x, ρ, u, v, p) = generateQuantities()
    writecsvstep(0, x, ρ, u, v, p)
    println("Step $(0)")
    gr()
    plot(x, ρ)

    # SPH operands
    L = zeros(PARTICLES_COUNT)
    H = zeros(PARTICLES_COUNT)

    # Simulation steps
    for step in 1:TIME_STEPS
        (x, ρ, u, v, p, L, H) = computeSimulationStep(step, x, ρ, u, v, p, L, H)
        writecsvstep(step, x, ρ, u, v, p)
        println("Step $(step)")
    end

    (x, ρ, u, v, p)
end
