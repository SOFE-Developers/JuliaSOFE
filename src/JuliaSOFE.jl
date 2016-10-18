__precompile__()

module JuliaSOFE

include(joinpath(dirname(@__FILE__), "elements", "Elements.jl"))
include(joinpath(dirname(@__FILE__), "meshes", "Meshes.jl"))

end # of module JuliaSOFE

