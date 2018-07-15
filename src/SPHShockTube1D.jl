module SPHShockTube1D

using Plots

include("kernels.jl")
include("output.jl")
include("initialconditions.jl")
include("solver.jl")
include("simulationstep.jl")

(x, ρ, u, v, p) = runsimulation()

gr()
plot(x, ρ)

directfind(x)

end # module
