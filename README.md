# EternityNum.lua

EternityNum.lua is a Lua library for representing extremely large numbers using layered exponential notation.

It supports values far beyond normal floating point limits, making it useful for:
- Incremental games
- Scientific simulations
- Hyper-operations
- Large-number mathematics

## Internal Representation

An EternityNum is stored as:

EN.new(sign, layer, mag)

Representing:

sign * 10^(10^(10^(...mag)))

where the exponent tower has `layer` levels.

Examples:
- Layer 0: `sign * mag`
- Layer 1: `sign * 10^mag`
- Layer 2: `sign * 10^(10^mag)`

Negative magnitudes at higher layers represent negative exponents.
