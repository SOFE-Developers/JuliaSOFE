__precompile__()

module JuliaSOFE

using Reexport

# Helpers
include(joinpath(dirname(@__FILE__), "helpers.jl"))
@reexport using .Helpers

# Elements
include(joinpath(dirname(@__FILE__), "elements", "Elements.jl"))
@reexport using .Elements

# Meshes
include(joinpath(dirname(@__FILE__), "meshes", "Meshes.jl"))
@reexport using .Meshes

# Quadrature
include(joinpath(dirname(@__FILE__), "quadrature", "QuadRule.jl"))
@reexport using .Quadrature

# Spaces
include(joinpath(dirname(@__FILE__), "spaces", "FESpace.jl"))
@reexport using .Spaces

# Operators
include(joinpath(dirname(@__FILE__), "operators", "Operators.jl"))
@reexport using .Operators

# Problems
include(joinpath(dirname(@__FILE__), "problems", "Problems.jl"))
@reexport using .Problems

# Preprocessing
include(joinpath(dirname(@__FILE__), "preprocessing", "MeshGeneration.jl"))
@reexport using .MeshGeneration

# Postprocessing
include(joinpath(dirname(@__FILE__), "postprocessing", "Visualization.jl"))
@reexport using .Visualization

end # of module JuliaSOFE

