module ThreeBodyDecaysIO

using ThreeBodyDecays
using ThreeBodyDecays: RecouplingLS, ParityRecoupling, NoRecoupling
using ThreeBodyDecays: AbstractDecayChain
using JSON
using OrderedCollections
using Parameters
using DataFrames

export serializeToDict
export topology2k
include("writer.jl")

export dict2kinematics
export dict2model
export dict2chain
# export dict2lineshape
# export dict2recoupling
include("reader.jl")

export BW
include("lineshapes.jl")

export string2complex
export update2values
include("utils.jl")

export validation_section
include("validation.jl")

end # module ThreeBodyDecaysIO
