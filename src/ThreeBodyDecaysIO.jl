module ThreeBodyDecaysIO

using ThreeBodyDecays
using ThreeBodyDecays: RecouplingLS, ParityRecoupling, NoRecoupling
using ThreeBodyDecays: AbstractDecayChain
using JSON
using OrderedCollections
using Parameters
using DataFrames
using HadronicLineshapes

export serializeToDict
export add_hs3_fields
include("writer.jl")

include("HadronicLineshapesIO.jl")

export dict2instance
# export dict2recoupling
include("reader.jl")

export BW
include("lineshapes.jl")

export reorder
export array2dict
export string2complex
export update2values
export topology2k, flatten_topology
include("utils.jl")

export validation_section
export angles_invariants
include("validation.jl")

end # module ThreeBodyDecaysIO
