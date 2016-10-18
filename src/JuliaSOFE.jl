__precompile__()

module JuliaSOFE

include(joinpath(dirname(@__FILE__), "elements", "Elements.jl"))
include(joinpath(dirname(@__FILE__), "meshes", "Meshes.jl"))
include(joinpath(dirname(@__FILE__), "quadrature", "QuadRule.jl"))

end # of module JuliaSOFE

