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

export HS3InputWrapper
include("HadronicLineshapesIO.jl")

export dict2kinematics
export dict2model
export dict2chain
export dict2lineshape
# export dict2recoupling
include("reader.jl")

export BW
include("lineshapes.jl")

export string2complex
export update2values
export topology2k, flatten_topology
include("utils.jl")

export validation_section
include("validation.jl")

end # module ThreeBodyDecaysIO
