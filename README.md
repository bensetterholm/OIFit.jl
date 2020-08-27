# OIFit.jl

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

A simple library for fitting geometrical models to OIFITS data, with routines originally implemented for [Setterholm, et al 2018](https://doi.org/10.3847/1538-4357/aaef2c).

## Installation

Installation is as simple as starting a Julia session, switching to the package manager mode (press the `]` key), and then entering:

```
add https://github.com/bensetterholm/OIFit.jl
```

This will install the OIFit.jl library and all dependencies. Once installed, the library can be updated anytime by entering
```
update
```
in package manager mode.

## Usage

```julia
using OIFit

OIFit.loaddir("path/to/dir") # Loads all .[oi]fits files into Vis2 structures
```

## License

[MIT](LICENSE)
