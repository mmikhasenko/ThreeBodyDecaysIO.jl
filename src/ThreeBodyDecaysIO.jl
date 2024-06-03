module ThreeBodyDecaysIO

using ThreeBodyDecays
using ThreeBodyDecays: RecouplingLS, ParityRecoupling, NoRecoupling
using ThreeBodyDecays: AbstractDecayChain
using JSON
using OrderedCollections
using Parameters
using DataFrames
using HadronicLineshapes
using MacroTools

export serializeToDict
export add_hs3_fields
export trivial_lineshape_parser
include("writer.jl")

export NamedArgFunc
include("HadronicLineshapesIO.jl")
export generic_function
include("generic_function.jl")

export HadronicUnpolarizedIntensity
include("HadronicUnpolarizedIntensity.jl")

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
export label_diff
include("utils.jl")

export angles_invariants
export validation_fields
export validation_section
include("validation.jl")

end # module ThreeBodyDecaysIO
