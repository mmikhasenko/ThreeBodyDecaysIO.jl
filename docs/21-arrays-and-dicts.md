# Arrays and Dictionaries: Design Choices in ThreeBodyDecaysIO

## JSON Arrays vs Julia Dictionaries

In serialization project, we made a deliberate design choice to use arrays for collection of objects,
to reassure a reproducible order of objects. For in-memory operations, we use dictionaries for fast lookups.

The `array2dict` function converts between these representations.
It is used in several JSON blocks, including:
- `variables` - for kinematic variable extraction (key: "node")
- `parameter_points` - for parameter point lookups (key: "name")
- `parameters` - for parameter value extraction (key: "name")
- `distributions` - for model content processing
