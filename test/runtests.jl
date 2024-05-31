using Test

# utils
include("test_utils.jl")

# lineshape reader
include("test_lineshape_reader.jl")
include("test_generic.jl")

# write_read_model
include("test_write_read_model.jl")

# test examples
include("test_model_content.jl")

# test validation
include("test_validation.jl")
