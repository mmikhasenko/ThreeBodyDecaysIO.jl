module ThreeBodyDecaysIO

using ThreeBodyDecays
using ThreeBodyDecays: RecouplingLS, ParityRecoupling, NoRecoupling
using ThreeBodyDecays: AbstractDecayChain
using JSON
using OrderedCollections
using Parameters


export wrap2dict
export valudation_section
include("dictwriters.jl")

end # module ThreeBodyDecaysIO
