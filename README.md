# ThreeBodyDecaysIO

## Overview

`ThreeBodyDecaysIO` is a Julia module designed for the serialization and deserialization of models related to three-body decays in particle physics. It leverages the `ThreeBodyDecays` framework to provide comprehensive support for reading and writing model descriptions, facilitating the exchange of complex decay models and their kinematic configurations.

## Features

-   Serialization and deserialization of three-body decay models to and from JSON format.
-   Support for detailed kinematics, lineshapes, and interaction chains descriptions.
-   Utilities for validation and manipulation of model components.

## Installation

To install `ThreeBodyDecaysIO`, use the Julia package manager. From the Julia REPL, type the following:

```julia
Pkg.add(url="https://github.com/mmikhasenko/ThreeBodyDecaysIO.jl")
```

## Usage

### Basic Example

Here is a basic example of using `ThreeBodyDecaysIO` to serialize a three-body decay model to JSON:

```julia
using ThreeBodyDecaysIO

# Define your model here
model = defineYourModel()

# Serialize to dictionary
dict = serializeToDict(model)

# Write to a JSON file
open("model.json", "w") do io
    JSON.print(io, dict, 4)
end
```

### Reading a Model

To read a model from a JSON file and parse it:

```julia
using ThreeBodyDecaysIO

# Read from JSON
json_content = open("model.json") do io
    JSON.parse(io)
end

# Parse the model
model = dict2instance(ThreeBodyDecay, json_content)
```

## Running Tests

To ensure `ThreeBodyDecaysIO` is working correctly, you can run its test suite:

```julia
using Pkg
Pkg.test("ThreeBodyDecaysIO")
```

This will execute a series of tests, verifying the functionality of model serialization/deserialization, kinematics parsing, and more. See details in [`test/runtests.jl`](test/runtests.jl)

## Running pre-commit hooks

To check code formatting, spelling, and other style issues before committing, you can run the pre-commit hooks manually:

```sh
pre-commit run --all-files
```

This will run all configured checks on the entire codebase and automatically fix some issues. See `.pre-commit-config.yaml` for details.

## Contributing

Contributions to `ThreeBodyDecaysIO` are welcome. To contribute, please fork the repository, make your changes, and submit a pull request.
