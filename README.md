# ParXRec.jl

[![build status](https://github.com/kamesy/ParXRec.jl/workflows/CI/badge.svg)](https://github.com/kamesy/ParXRec.jl/actions?query=workflow%3ACI)

Julia package for loading Philips' PAR/REC and XML/REC files.

## Installation

ParXRec requires Julia v1.0 or later.

```julia
julia> ]add ParXRec
```

## Usage

Load PAR/XML header and REC data:
```julia
using ParXRec

rec = ParXRec.load("file.{PAR, XML, REC}")
hdr = rec.hdr
data = rec.data
```

Load PAR/XML header only:
```julia
rec = ParXRec.load("file.{PAR, XML, REC}", load_data = false)
hdr = rec.hdr
```

Extract voxel size and echo times from header:
```julia
vsz = voxelsize(hdr)
TEs = echotimes(hdr)
```

## Notes
- Supports PAR versions 3, 4, 4.1, 4.2.
- REC data is sorted and loaded into a single array. The dimensions are:
    - slices
    - echoes
    - dynamics
    - cardiac phases
    - diffusion b values
    - gradient orientations
    - label types
    - image types
    - sequences

    Singleton dimensions (except for slices) are removed.
