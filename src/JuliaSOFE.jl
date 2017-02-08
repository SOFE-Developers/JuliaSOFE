__precompile__()

module JuliaSOFE

using Reexport

# Helpers
include(joinpath(dirname(@__FILE__), "helpers.jl"))
@reexport using .Helpers

# Meshes
# include(joinpath(dirname(@__FILE__), "meshes", "Meshes.jl"))
# @reexport using .Meshes

# CMesh Backend
@reexport using CMesh

# Elements
include(joinpath(dirname(@__FILE__), "elements", "Elements.jl"))
@reexport using .Elements

# CMeshExtensions
include(joinpath(dirname(@__FILE__), "cmeshext", "CMeshExtensions.jl"))
@reexport using .CMeshExtensions

# Quadrature
include(joinpath(dirname(@__FILE__), "quadrature", "QuadRule.jl"))
@reexport using .Quadrature

# Spaces
include(joinpath(dirname(@__FILE__), "spaces", "FESpace.jl"))
@reexport using .Spaces

# Operators
include(joinpath(dirname(@__FILE__), "operators", "Operators.jl"))
@reexport using .Operators

# Problems / Constraints
include(joinpath(dirname(@__FILE__), "problems", "Problems.jl"))
include(joinpath(dirname(@__FILE__), "problems", "Constraints.jl"))
@reexport using .Problems
@reexport using .Constraints

# Postprocessing
include(joinpath(dirname(@__FILE__), "postprocessing", "Visualization.jl"))
@reexport using .Visualization

end # of module JuliaSOFE

