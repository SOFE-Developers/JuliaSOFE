__precompile__()

module JuliaSOFE

include(joinpath(dirname(@__FILE__), "helpers.jl"))
using .Helpers

include(joinpath(dirname(@__FILE__), "elements", "Elements.jl"))
using .Elements
export LagrangeP1, LagrangeQ1

include(joinpath(dirname(@__FILE__), "meshes", "Meshes.jl"))
using .Meshes
export Mesh, TensorProductMesh, UnitSquare, UnitCube
export getNodes, getConnectivity, getEntities, getNumber

include(joinpath(dirname(@__FILE__), "quadrature", "QuadRule.jl"))
using .Quadrature
export QuadRule, QuadRuleSimp1

include(joinpath(dirname(@__FILE__), "spaces", "FESpace.jl"))
using .Spaces
export FESpace, dofMap

include(joinpath(dirname(@__FILE__), "operators", "Operators.jl"))
using .Operators

end # of module JuliaSOFE

