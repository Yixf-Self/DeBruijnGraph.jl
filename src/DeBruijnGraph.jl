module DeBruijnGraph

# package code goes here

using Graphs

export Node, DebruijnGraph,

       isSemiBalanced, isBalanced, name,

       nnodes,nedges,hasEulerianWalk,hasEulerianCycle,isEulerian,eulerianWalkOrCyle




include("graph/graph.jl")
include("visual/visual.jl")
include("assembly/assembly.jl")

end # module
