using SPHShockTube1D
# @static if VERSION < v"0.7.0-DEV.2005"
#     using Base.Test
# else
#     using Test
# end

# write your own tests here
(x, ρ, u, v, p) = SPHShockTube1D.runsimulation()
gr()
plot(x, ρ)
