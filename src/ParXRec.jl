module ParXRec


using AxisArrays
using Dates
using LightXML
using Printf


export readrec, readpar, readxml, echotimes, voxelsize, recdims


###
### Types
###

const UNKNOWN_KEY = :Other


"""
    mutable struct SeriesInfo
        Patient_Name            :: String            = ""
        Examination_Name        :: String            = ""
        Protocol_Name           :: String            = ""
        Examination_Date        :: Date              = Date(1970, 1, 1)
        Examination_Time        :: Time              = Time(0, 0, 0)
        Series_Data_Type        :: String            = ""
        Aquisition_Number       :: Int               = 0
        Reconstruction_Number   :: Int               = 0
        Scan_Duration           :: Float64           = 0
        Max_No_Phases           :: Int               = 1
        Max_No_Echoes           :: Int               = 1
        Max_No_Slices           :: Int               = 1
        Max_No_Dynamics         :: Int               = 1
        Max_No_Mixes            :: Int               = 1
        Max_No_B_Values         :: Int               = 1
        Max_No_Gradient_Orients :: Int               = 1
        No_Label_Types          :: Int               = 0
        Patient_Position        :: String            = ""
        Preparation_Direction   :: String            = ""
        Technique               :: String            = ""
        Scan_Resolution_X       :: Int               = 0
        Scan_Resolution_Y       :: Int               = 0
        Scan_Mode               :: String            = ""
        Repetition_Times        :: NTuple{2,Float64} = (0, 0)
        FOV_AP                  :: Float64           = 0
        FOV_FH                  :: Float64           = 0
        FOV_RL                  :: Float64           = 0
        Water_Fat_Shift         :: Float64           = 0
        Angulation_AP           :: Float64           = 0
        Angulation_FH           :: Float64           = 0
        Angulation_RL           :: Float64           = 0
        Off_Center_AP           :: Float64           = 0
        Off_Center_FH           :: Float64           = 0
        Off_Center_RL           :: Float64           = 0
        Flow_Compensation       :: Bool              = false
        Presaturation           :: Bool              = false
        Phase_Encoding_Velocity :: NTuple{3,Float64} = (0, 0, 0)
        MTC                     :: Bool              = false
        SPIR                    :: Bool              = false
        EPI_factor              :: Int               = 0
        Dynamic_Scan            :: Bool              = false
        Diffusion               :: Bool              = false
        Diffusion_Echo_Time     :: Float64           = 0
        Other                   :: Dict{String,Any}  = Dict()
    end

Struct containing .PAR or .XML series info.
"""
Base.@kwdef mutable struct SeriesInfo
    Patient_Name            :: String            = ""
    Examination_Name        :: String            = ""
    Protocol_Name           :: String            = ""
    Examination_Date        :: Date              = Date(1970, 1, 1)
    Examination_Time        :: Time              = Time(0, 0, 0)
    Series_Data_Type        :: String            = ""
    Aquisition_Number       :: Int               = 0
    Reconstruction_Number   :: Int               = 0
    Scan_Duration           :: Float64           = 0
    Max_No_Phases           :: Int               = 1
    Max_No_Echoes           :: Int               = 1
    Max_No_Slices           :: Int               = 1
    Max_No_Dynamics         :: Int               = 1
    Max_No_Mixes            :: Int               = 1
    Max_No_B_Values         :: Int               = 1
    Max_No_Gradient_Orients :: Int               = 1
    No_Label_Types          :: Int               = 0
    Patient_Position        :: String            = ""
    Preparation_Direction   :: String            = ""
    Technique               :: String            = ""
    Scan_Resolution_X       :: Int               = 0
    Scan_Resolution_Y       :: Int               = 0
    Scan_Mode               :: String            = ""
    Repetition_Times        :: NTuple{2,Float64} = (0, 0)
    FOV_AP                  :: Float64           = 0
    FOV_FH                  :: Float64           = 0
    FOV_RL                  :: Float64           = 0
    Water_Fat_Shift         :: Float64           = 0
    Angulation_AP           :: Float64           = 0
    Angulation_FH           :: Float64           = 0
    Angulation_RL           :: Float64           = 0
    Off_Center_AP           :: Float64           = 0
    Off_Center_FH           :: Float64           = 0
    Off_Center_RL           :: Float64           = 0
    Flow_Compensation       :: Bool              = false
    Presaturation           :: Bool              = false
    Phase_Encoding_Velocity :: NTuple{3,Float64} = (0, 0, 0)
    MTC                     :: Bool              = false
    SPIR                    :: Bool              = false
    EPI_factor              :: Int               = 0
    Dynamic_Scan            :: Bool              = false
    Diffusion               :: Bool              = false
    Diffusion_Echo_Time     :: Float64           = 0
    Other                   :: Dict{String,Any}  = Dict()
end


"""
    mutable struct ImageKey
        Slice       :: Int    = 1
        Echo        :: Int    = 1
        Dynamic     :: Int    = 1
        Phase       :: Int    = 1
        BValue      :: Int    = 1
        Grad_Orient :: Int    = 1
        Label_Type  :: String = "" # Enumeration
        Type        :: String = "" # Enumeration
        Sequence    :: String = "" # Enumeration
        Index       :: Int    = 0
    end

Struct containing .PAR or .XML image keys.
"""
Base.@kwdef mutable struct ImageKey
    Slice       :: Int    = 1
    Echo        :: Int    = 1
    Dynamic     :: Int    = 1
    Phase       :: Int    = 1
    BValue      :: Int    = 1
    Grad_Orient :: Int    = 1
    Label_Type  :: String = "" # Enumeration
    Type        :: String = "" # Enumeration
    Sequence    :: String = "" # Enumeration
    Index       :: Int    = 0
end


"""
    mutable struct ImageInfo
        Key                        :: ImageKey
        Pixel_Size                 :: Int               = 0
        Scan_Percentage            :: Float64           = 0
        Resolution_X               :: Int               = 0
        Resolution_Y               :: Int               = 0
        Rescale_Intercept          :: Float64           = 0
        Rescale_Slope              :: Float64           = 0
        Scale_Slope                :: Float64           = 0
        Window_Center              :: Float64           = 0
        Window_Width               :: Float64           = 0
        Slice_Thickness            :: Float64           = 0
        Slice_Gap                  :: Float64           = 0
        Display_Orientation        :: String            = ""      # Enumeration
        fMRI_Status_Indication     :: Int               = 0
        Image_Type_Ed_Es           :: String            = ""      # Enumeration
        Pixel_Spacing              :: NTuple{2,Float64} = (0, 0)
        Echo_Time                  :: Float64           = 0
        Dyn_Scan_Begin_Time        :: Float64           = 0
        Trigger_Time               :: Float64           = 0
        Diffusion_B_Factor         :: Float64           = 0
        No_Averages                :: Float64           = 0
        Image_Flip_Angle           :: Float64           = 0
        Cardiac_Frequency          :: Int               = 0
        Min_RR_Interval            :: Int               = 0
        Max_RR_Interval            :: Int               = 0
        TURBO_Factor               :: Int               = 0
        Inversion_Delay            :: Float64           = 0
        Contrast_Type              :: String            = ""
        Diffusion_Anisotropy_Type  :: String            = ""
        Diffusion_AP               :: Float64           = 0
        Diffusion_FH               :: Float64           = 0
        Diffusion_RL               :: Float64           = 0
        Angulation_AP              :: Float64           = 0
        Angulation_FH              :: Float64           = 0
        Angulation_RL              :: Float64           = 0
        Offcenter_AP               :: Float64           = 0
        Offcenter_FH               :: Float64           = 0
        Offcenter_RL               :: Float64           = 0
        Slice_Orientation          :: String            = ""      # Enumeration
        Image_Planar_Configuration :: Int               = 0
        Samples_Per_Pixel          :: Int               = 0
        Other                      :: Dict{String,Any}  = Dict()
    end

Struct containing .PAR or .XML image info.
"""
Base.@kwdef mutable struct ImageInfo
    Key                        :: ImageKey
    Pixel_Size                 :: Int               = 0
    Scan_Percentage            :: Float64           = 0
    Resolution_X               :: Int               = 0
    Resolution_Y               :: Int               = 0
    Rescale_Intercept          :: Float64           = 0
    Rescale_Slope              :: Float64           = 0
    Scale_Slope                :: Float64           = 0
    Window_Center              :: Float64           = 0
    Window_Width               :: Float64           = 0
    Slice_Thickness            :: Float64           = 0
    Slice_Gap                  :: Float64           = 0
    Display_Orientation        :: String            = ""      # Enumeration
    fMRI_Status_Indication     :: Int               = 0
    Image_Type_Ed_Es           :: String            = ""      # Enumeration
    Pixel_Spacing              :: NTuple{2,Float64} = (0, 0)
    Echo_Time                  :: Float64           = 0
    Dyn_Scan_Begin_Time        :: Float64           = 0
    Trigger_Time               :: Float64           = 0
    Diffusion_B_Factor         :: Float64           = 0
    No_Averages                :: Float64           = 0
    Image_Flip_Angle           :: Float64           = 0
    Cardiac_Frequency          :: Int               = 0
    Min_RR_Interval            :: Int               = 0
    Max_RR_Interval            :: Int               = 0
    TURBO_Factor               :: Int               = 0
    Inversion_Delay            :: Float64           = 0
    Contrast_Type              :: String            = ""
    Diffusion_Anisotropy_Type  :: String            = ""
    Diffusion_AP               :: Float64           = 0
    Diffusion_FH               :: Float64           = 0
    Diffusion_RL               :: Float64           = 0
    Angulation_AP              :: Float64           = 0
    Angulation_FH              :: Float64           = 0
    Angulation_RL              :: Float64           = 0
    Offcenter_AP               :: Float64           = 0
    Offcenter_FH               :: Float64           = 0
    Offcenter_RL               :: Float64           = 0
    Slice_Orientation          :: String            = ""      # Enumeration
    Image_Planar_Configuration :: Int               = 0
    Samples_Per_Pixel          :: Int               = 0
    Other                      :: Dict{String,Any}  = Dict()
end


"""
    struct RecHeader
        series :: SeriesInfo
        images :: Vector{ImageInfo}
    end

Struct containing .PAR or .XML header.
"""
struct RecHeader
    series :: SeriesInfo
    images :: Vector{ImageInfo}
end


"""
    struct Rec{A<:Union{Nothing, AbstractArray}}
        hdr  :: RecHeader
        data :: A
    end

Struct containing .PAR or .XML header and .REC data array (or nothing).
"""
struct Rec{A<:Union{Nothing, AbstractArray}}
    hdr  :: RecHeader
    data :: A
end


###
### Interface
###

"""
    load(path; kwargs...) = load(Float32, path; kwargs...)

    load(
        ::Type{T},
        path::AbstractString;
        load_data::Bool = true,
        scale_data::Bool = true
    ) where {T<:AbstractFloat} -> Rec{Union{Nothing, AxisArray}}

Load Philips PAR/REC or XML/REC header and data array.

### Arguments
- `path::AbstractString`: path to .REC or .XML or .PAR file

### Keywords
- `load_data::Bool = true`: load data and header (true) or header only (false)
- `scale_data::Bool = true`: scale data array with
    `(data + rescale slope * rescale intercept) / scale slope` (true) or
    return raw array (false)

### Returns
- `Rec{Union{Nothing, AxisArray}}`: struct containing header and data array
    if `load_data = true` or nothing if `load_data = false`
"""
load(path; kwargs...) = load(Float32, path; kwargs...)

function load(
    ::Type{T},
    path::AbstractString;
    load_data::Bool = true,
    scale_data::Bool = true,
) where {T<:AbstractFloat}
    hdrpath, recpath = parxrec_paths(path, load_data)
    _, ext = splitext(hdrpath)

    hdr = lowercase(ext) == ".par" ?
        readpar(hdrpath) :
        readxml(hdrpath)

    data = load_data ?
        readrec(T, recpath, hdr, scale_data=scale_data) :
        nothing

    return Rec(hdr, data)
end


###
### .REC
###

const _RECPREC = Dict{Int, DataType}(
    8  => UInt8,
    16 => UInt16,
    32 => Float32,
    64 => Float64,
)

"""
    readrec(path, hdr; kwargs...) = readrec(Float32, path, hdr; kwargs...)

    function readrec(
        ::Type{T},
        path::AbstractString,
        hdr::RecHeader;
        scale_data::Bool = true
    ) where {T<:AbstractFloat} -> AxisArray

Load .REC file.

### Arguments
- `path::AbstractString`: path to .REC file
- `hdr::RecHeader`: struct containing .PAR or .XML header

### Keywords
- `scale_data::Bool = true`: scale data array with
    `(data + rescale slope * rescale intercept) / scale slope` (true) or
    return raw array (false)

### Returns
- `AxisArray`: .REC data array
"""
readrec(path, hdr; kwargs...) = readrec(Float32, path, hdr; kwargs...)

function readrec(
    ::Type{T},
    path::AbstractString,
    hdr::RecHeader;
    scale_data::Bool = true
) where {T<:AbstractFloat}
    validate(hdr)

    imgs = sort(hdr.images, by = (i -> i.Key.Index))

    ps = Int(imgs[1].Pixel_Size)
    TR = _RECPREC[ps]

    dims = recdims(hdr)

    ax0 = (
        nx       = maximum(dims[:nx]),
        ny       = maximum(dims[:ny]),
        nz       = maximum(dims[:nz]),
        echoes   = maximum(dims[:echoes]),
        dynamics = maximum(dims[:dynamics]),
        phases   = maximum(dims[:phases]),
        bvals    = maximum(dims[:bvals]),
        grads    = maximum(dims[:grads]),
        labels   = length(dims[:labels]),
        types    = length(dims[:types]),
        seqs     = length(dims[:seqs]),
    )

    # remove singleton dimensions but keep nx, ny, nz
    ax = (; [(a, n) for (a, n) in pairs(ax0) if n > 1 || a ∈ (:nx, :ny, :nz)]...)

    nx = ax[:nx]
    ny = ax[:ny]

    labels = Dict(k => i for (i, k) in enumerate(dims[:labels]))
    types = Dict(k => i for (i, k) in enumerate(dims[:types]))
    seqs = Dict(k => i for (i, k) in enumerate(dims[:seqs]))

    L = LinearIndices(values(ax0))
    data = zeros(T, values(ax))
    slice = zeros(TR, nx * ny)
    slicet = zeros(TR, nx, ny)

    nb = nx * ny * ps÷8

    rsz = filesize(path)
    w = (nb * length(imgs)) - rsz
    w < 0 && @warn "REC file larger than expected from hdr"
    w > 0 && @warn "REC file smaller than expected from hdr"
    rsz -= nb

    open(path, "r") do io
        @inbounds for i in imgs
            pos = i.Key.Index * nb
            pos > rsz && continue

            seek(io, pos)
            read!(io, slicet)
            copyto!(slice, transpose(slicet))

            i3  = i.Key.Slice
            i4  = i.Key.Echo
            i5  = i.Key.Dynamic
            i6  = i.Key.Phase
            i7  = i.Key.BValue
            i8  = i.Key.Grad_Orient
            i9  = labels[i.Key.Label_Type]
            i10 = types[(i.Key.Type, i.Rescale_Intercept)]
            i11 = seqs[i.Key.Sequence]

            I = L[1,1,i3,i4,i5,i6,i7,i8,i9,i10,i11] : L[nx,ny,i3,i4,i5,i6,i7,i8,i9,i10,i11]

            if scale_data
                ss = i.Scale_Slope ≈ 0.0 ? 1.0 : inv(i.Scale_Slope)
                rs = i.Rescale_Slope ≈ 0.0 ? 1.0 : inv(i.Rescale_Slope)
                ri = rs * i.Rescale_Intercept
                data[I] .= T(ss) .* (slice .+ T(ri))
            else
                data[I] .= slice
            end
        end
    end

    return AxisArray(data, keys(ax)...)
end


###
### .XML
###

const _XML_HEADER           = "PRIDE_V5"
const _XML_SERIES_HEADER    = "Series_Info"
const _XML_IMAGE_ARR_HEADER = "Image_Array"
const _XML_IMAGE_HEADER     = "Image_Info"
const _XML_IMAGE_KEY_HEADER = "Key"
const _XML_ATTRIB_HEADER    = "Attribute"


"""
    readxml(path::AbstractString) -> RecHeader

Load .XML header.

### Arguments
- `path::AbstractString`: path to .XML file

### Returns
- `RecHeader`: struct containing .XML header
"""
function readxml(path::AbstractString)
    series = SeriesInfo()
    images = Vector{ImageInfo}()

    Ts = NamedTuple{fieldnames(SeriesInfo)}(SeriesInfo.types)
    Tk = NamedTuple{fieldnames(ImageKey)}(ImageKey.types)
    Ti = NamedTuple{fieldnames(ImageInfo)}(ImageInfo.types)

    xdoc = parse_file(path)
    xroot = root(xdoc)
    name(xroot) != _XML_HEADER &&
        error("Unknown XML/REC header: $(name(xroot)) (expected <$(_XML_HEADER)>)")

    # Series Info
    cur = find_element(xroot, _XML_SERIES_HEADER)
    cur === nothing &&
        error("<$(_XML_SERIES_HEADER)> tag not found")

    for e in child_elements(cur)
        name(e) != _XML_ATTRIB_HEADER && continue
        sethdr!(series, Ts, e)
    end

    # Image Info
    cur = find_element(xroot, _XML_IMAGE_ARR_HEADER)
    cur === nothing &&
        error("<$(_XML_IMAGE_ARR_HEADER)> tag not found")

    for c in child_elements(cur)
        name(c) != _XML_IMAGE_HEADER && continue

        ks = find_element(c, _XML_IMAGE_KEY_HEADER)
        ks === nothing &&
            error("<$(_XML_IMAGE_KEY_HEADER)> tag not found")

        imagekey = ImageKey()
        for e in child_elements(ks)
            name(e) != _XML_ATTRIB_HEADER && continue
            sethdr!(imagekey, Tk, e)
        end

        image = ImageInfo(Key = imagekey)
        for e in child_elements(c)
            name(e) != _XML_ATTRIB_HEADER && continue
            sethdr!(image, Ti, e)
        end

        push!(images, image)
    end

    free(xdoc)
    return RecHeader(series, images)
end


function sethdr!(h, T::NamedTuple, e::XMLElement)
    name = attribute(e, "Name")
    name === nothing && return nothing

    k = Symbol(replace(name, " " => "_"))
    v = content(e)

    if !haskey(T, k)
        _u = getfield(h, UNKNOWN_KEY)
        _u[name] = v
    elseif v != ""
        setfield!(h, k, _parse(T[k], v))
    end

    return h
end


###
### .PAR
###

const _PARSERIES = Dict{String, Union{Symbol, Tuple}}(
    "Patient name"                        => :Patient_Name,
    "Examination name"                    => :Examination_Name,
    "Protocol name"                       => :Protocol_Name,
    "Examination date/time"               => (:Examination_Date, :Examination_Time),
    "Series Type"                         => :Series_Data_Type,
    "Acquisition nr"                      => :Aquisition_Number,
    "Reconstruction nr"                   => :Reconstruction_Number,
    "Scan Duration [sec]"                 => :Scan_Duration,
    "Max. number of cardiac phases"       => :Max_No_Phases,
    "Max. number of echoes"               => :Max_No_Echoes,
    "Max. number of slices/locations"     => :Max_No_Slices,
    "Max. number of dynamics"             => :Max_No_Dynamics,
    "Max. number of mixes"                => :Max_No_Mixes,
    "Patient position"                    => :Patient_Position,
    "Preparation direction"               => :Preparation_Direction,
    "Technique"                           => :Technique,
    "Scan resolution  (x, y)"             => (:Scan_Resolution_X, :Scan_Resolution_Y),
    "Scan mode"                           => :Scan_Mode,
    "Repetition time [ms]"                => :Repetition_Times,
    "FOV (ap,fh,rl) [mm]"                 => (:FOV_AP, :FOV_FH, :FOV_RL),
    "Water Fat shift [pixels]"            => :Water_Fat_Shift,
    "Angulation midslice(ap,fh,rl)[degr]" => (:Angulation_AP, :Angulation_FH, :Angulation_RL),
    "Off Centre midslice(ap,fh,rl) [mm]"  => (:Off_Center_AP, :Off_Center_FH, :Off_Center_RL),
    "Flow compensation <0=no 1=yes> ?"    => :Flow_Compensation,
    "Presaturation     <0=no 1=yes> ?"    => :Presaturation,
    "Phase encoding velocity [cm/sec]"    => :Phase_Encoding_Velocity,
    "MTC               <0=no 1=yes> ?"    => :MTC,
    "SPIR              <0=no 1=yes> ?"    => :SPIR,
    "EPI factor        <0,1=no EPI>"      => :EPI_factor,
    "Dynamic scan      <0=no 1=yes> ?"    => :Dynamic_Scan,
    "Diffusion         <0=no 1=yes> ?"    => :Diffusion,
    "Diffusion echo time [ms]"            => :Diffusion_Echo_Time,
    "Max. number of diffusion values"     => :Max_No_B_Values,
    "Max. number of gradient orients"     => :Max_No_Gradient_Orients,
    "Number of label types   <0=no ASL>"  => :No_Label_Types,
)


const _PARSERIES3 = Dict{String, Union{Symbol, Tuple}}(
    "Image pixel size [8 or 16 bits]"     => :Pixel_Size,
    "Scan percentage"                     => :Scan_Percentage,
    "Recon resolution (x, y)"             => (:Resolution_X, :Resolution_Y),
    "Slice thickness [mm]"                => :Slice_Thickness,
    "Slice gap [mm]"                      => :Slice_Gap,
    "Number of averages"                  => :No_Averages,
    "Cardiac frequency"                   => :Cardiac_Frequency,
    "Min. RR interval"                    => :Min_RR_Interval,
    "Max. RR interval"                    => :Max_RR_Interval,
    "TURBO factor      <0=no turbo>"      => :TURBO_Factor,
    "Inversion delay [ms]"                => :Inversion_Delay,
)


const _PARKEY = Dict{Int, Symbol}(
    1  => :Slice,
    2  => :Echo,
    3  => :Dynamic,
    4  => :Phase,
    5  => :Type,
    6  => :Sequence,
    7  => :Index,
    42 => :BValue,
    43 => :Grad_Orient,
    49 => :Label_Type,
)


const _PARIMAGE = Dict{Int, Symbol}(
    8  => :Pixel_Size,
    9  => :Scan_Percentage,
    10 => :Resolution_X,
    11 => :Resolution_Y,
    12 => :Rescale_Intercept,
    13 => :Rescale_Slope,
    14 => :Scale_Slope,
    15 => :Window_Center,
    16 => :Window_Width,
    17 => :Angulation_AP,
    18 => :Angulation_FH,
    19 => :Angulation_RL,
    20 => :Offcenter_AP,
    21 => :Offcenter_FH,
    22 => :Offcenter_RL,
    23 => :Slice_Thickness,
    24 => :Slice_Gap,
    25 => :Display_Orientation,
    26 => :Slice_Orientation,
    27 => :fMRI_Status_Indication,
    28 => :Image_Type_Ed_Es,
    29 => :Pixel_Spacing,  # [29, 30]
    31 => :Echo_Time,
    32 => :Dyn_Scan_Begin_Time,
    33 => :Trigger_Time,
    34 => :Diffusion_B_Factor,
    35 => :No_Averages,
    36 => :Image_Flip_Angle,
    37 => :Cardiac_Frequency,
    38 => :Min_RR_Interval,
    39 => :Max_RR_Interval,
    40 => :TURBO_Factor,
    41 => :Inversion_Delay,
    44 => :Contrast_Type,
    45 => :Diffusion_Anisotropy_Type,
    46 => :Diffusion_AP,
    47 => :Diffusion_FH,
    48 => :Diffusion_RL,
)


const _PARIMAGE3 = Dict{Int, Symbol}(
    8  => :Rescale_Intercept,
    9  => :Rescale_Slope,
    10 => :Scale_Slope,
    11 => :Window_Center,
    12 => :Window_Width,
    13 => :Angulation_AP,
    14 => :Angulation_FH,
    15 => :Angulation_RL,
    16 => :Offcenter_AP,
    17 => :Offcenter_FH,
    18 => :Offcenter_RL,
    19 => :Display_Orientation,
    20 => :Slice_Orientation,
    21 => :fMRI_Status_Indication,
    22 => :Image_Type_Ed_Es,
    23 => :Pixel_Spacing,  # [23, 24]
    25 => :Echo_Time,
    26 => :Dyn_Scan_Begin_Time,
    27 => :Trigger_Time,
    28 => :Diffusion_B_Factor,
    29 => :Image_Flip_Angle,
)


"""
    readpar(path::AbstractString) -> RecHeader

Load .PAR header.

### Arguments
- `path::AbstractString`: path to .PAR file

### Returns
- `RecHeader`: struct containing .PAR header
"""
function readpar(path::AbstractString)
    version = ""
    gyrotools = false

    series = SeriesInfo()
    images = Vector{ImageInfo}()

    Ts = NamedTuple{fieldnames(SeriesInfo)}(SeriesInfo.types)
    Tk = NamedTuple{fieldnames(ImageKey)}(ImageKey.types)
    Ti = NamedTuple{fieldnames(ImageInfo)}(ImageInfo.types)

    doc = eachline(path)
    iter = iterate(doc)

    # PAR version number
    while iter !== nothing
        line, state = iter
        first(line) == '.' && break

        if occursin("CLINICAL TRYOUT", line)
            m = match(r"(?<=V)[0-9]\.?[0-9]*$", strip(line))
            if m !== nothing
                version = m.match
            end
            iter = iterate(doc, state)
            break

        elseif occursin(r"GyroTools|gyrotools", line)
            gyrotools = true
        end

        iter = iterate(doc, state)
    end

    v = tryparse(Float64, version)
    if v === nothing
        error("PAR version not found")

    elseif v < 3
        doc = nothing
        iter = nothing
        error("PAR V$(version) not supported")

    elseif v < 4
        v3 = true
        image3 = ImageInfo(Key = ImageKey())
        parimage = _PARIMAGE3

    else
        v3 = false
        parimage = _PARIMAGE
    end

    # Series Info / General Information
    while iter !== nothing
        line, state = iter

        if line == ""
            iter = iterate(doc, state)
            break
        end

        if line[1] != '.'
            (line[1] == ' ' || isdigit(line[1])) && break
            iter = iterate(doc, state)
            continue
        end

        if line[end] == ':'
            iter = iterate(doc, state)
            continue
        end

        ks, vs = strip.(split(line[2:end], ": "))

        if gyrotools && ks == "Patient Position"
            ks = "Patient position"
        end

        if (gyrotools || v3) && occursin("[msec]", ks)
            ks = replace(ks, "[msec]" => "[ms]")
        end

        if haskey(_PARSERIES, ks)
            k = _PARSERIES[ks]

            if k == :Repetition_Times
                if length(split(vs)) == 1
                    vs *= " 0.0"
                end
            end

            if k isa Symbol
                setfield!(series, k, _parse(Ts[k], vs))

            elseif k isa Tuple
                if occursin("Examination", ks)
                    x = split(vs, "/")
                    # in case month and/or day is 00
                    x[1] = replace(x[1], ".00" => ".01")
                else
                    x = split(vs)
                end

                x = strip.(x)

                for (i, _k) in enumerate(k)
                    setfield!(series, _k, _parse(Ts[_k], x[i]))
                end

            else
                error("Not implemented")
            end

        elseif v3 && haskey(_PARSERIES3, ks)
            k = _PARSERIES3[ks]

            if k isa Symbol
                setfield!(image3, k, _parse(Ti[k], vs))

            elseif k isa Tuple
                x = strip.(split(vs))
                for (i, _k) in enumerate(k)
                    setfield!(image3, _k, _parse(Ti[_k], x[i]))
                end

            else
                error("Not implemented")
            end

        else
            _u = getfield(series, UNKNOWN_KEY)
            _u[ks] = vs
        end

        iter = iterate(doc, state)
    end

    # Image Info
    while iter !== nothing
        line, state = iter

        if line == "" || (line[1] != ' ' && !isdigit(line[1]))
            iter = iterate(doc, state)
            continue
        end

        xs = split(line)
        N = length(xs)

        imagekey = ImageKey()
        for (i, k) in _PARKEY
            i > N && continue
            setfield!(imagekey, k, _parse(Tk[k], xs[i]))
        end

        image = ImageInfo(Key = imagekey)
        for (i, k) in parimage
            i > N && continue
            k == :Pixel_Spacing ?
                setfield!(image, k, _parse(Ti[k], xs[i]*" "*xs[i+1])) :
                setfield!(image, k, _parse(Ti[k], xs[i]))
        end

        if v3
            for v in values(_PARSERIES3)
                if v isa Tuple
                    for _v in v
                        setfield!(image, _v, getfield(image3, _v))
                    end
                else
                    setfield!(image, v, getfield(image3, v))
                end
            end
        end

        push!(images, image)
        iter = iterate(doc, state)
    end

    return RecHeader(series, images)
end


###
### Utility
###

"""
    echotimes(R::Rec)
    echotimes(H::RecHeader)
    echotimes(I::AbstractArray{ImageInfo}) -> NTuple{N, Float64}

Get echo times (in seconds).
"""
echotimes(R::Rec) = echotimes(R.hdr)
echotimes(H::RecHeader) = echotimes(H.images)
echotimes(I::AbstractArray{ImageInfo}) =
    1e-3 .* (sort!(unique!((i -> i.Echo_Time).(I)))...,)


"""
    voxelsize(R::Rec)
    voxelsize(H::RecHeader)
    voxelsize(I::ImageInfo) -> NTuple{3, Float64}
    voxelsize(I::AbstractArray{ImageInfo}) -> Union{NTuple{3, Float64}, Array{NTuple{3, Float64}}}

Get voxel size.
"""
voxelsize(R::Rec) = voxelsize(R.hdr)
voxelsize(H::RecHeader) = voxelsize(H.images)
voxelsize(I::ImageInfo) = (I.Pixel_Spacing..., I.Slice_Thickness + I.Slice_Gap)
voxelsize(I::AbstractArray{ImageInfo}) =
    unique!(voxelsize.(I)) |> v -> length(v) == 1 ? first(v) : v


"""
    recdims(H::RecHeader) -> Dict{Symbol, Int}

Array dimensions for REC data array from .PAR/.XML header. Dict is unordered.

### Keys
- `:nx`
- `:ny`
- `:nz`
- `:echoes`
- `:dynamics`
- `:phases`
- `:bvals`
- `:grads`
- `:labels`
- `:types`
- `:seqs`
"""
function recdims(H::RecHeader)
    K = (i -> i.Key).(H.images)

    @inline fun(f, k; by=identity) = sort!(unique!(f.(k)), by=by)

    nx       = fun(i -> i.Resolution_X, H.images)
    ny       = fun(i -> i.Resolution_Y, H.images)
    nz       = fun(k -> k.Slice, K)

    echoes   = fun(k -> k.Echo, K)
    dynamics = fun(k -> k.Dynamic, K)
    phases   = fun(k -> k.Phase, K)
    bvals    = fun(k -> k.BValue, K)
    grads    = fun(k -> k.Grad_Orient, K)

    labels   = fun(k -> k.Label_Type, K)
    seqs     = fun(k -> k.Sequence, K)

    # derived images have type = -1 in .PAR.
    # differentiate based on rescale intercept.
    types    = fun(
        i -> (i.Key.Type, i.Rescale_Intercept),
        H.images,
        by = t -> (t[1], abs(t[2])) # phase -3142, mag 0 -> flip order
    )

    return Dict(
        :nx       => nx,
        :ny       => ny,
        :nz       => nz,

        :echoes   => echoes,
        :dynamics => dynamics,
        :phases   => phases,
        :bvals    => bvals,
        :grads    => grads,

        :labels   => labels,
        :types    => types,
        :seqs     => seqs,
    )
end


###
### Misc
###

function parxrec_paths(path::AbstractString, load_data::Bool)
    !isfile(path) && error("File `$path` does not exist")

    path1, ext = splitext(path)
    ext = lowercase(ext)

    if ext == ".xml" || ext == ".par"
        hdrpath = path
        recpath = load_data ? (
            isfile(path1 * ".REC") ? path1 * ".REC" :
            isfile(path1 * ".rec") ? path1 * ".rec" :
            error(".REC file $path1 does not exist")
        ) : nothing

    elseif ext == ".rec"
        recpath = path
        hdrpath =
            isfile(path1 * ".XML") ? path1 * ".XML" :
            isfile(path1 * ".PAR") ? path1 * ".PAR" :
            isfile(path1 * ".xml") ? path1 * ".xml" :
            isfile(path1 * ".par") ? path1 * ".par" :
            error(".PAR/.XML file $path1 does not exist")
    else
        error("File $path does not have valid PAR/XML/REC extension")
    end

    return hdrpath, recpath
end


function validate(hdr::RecHeader)
    series = hdr.series
    images = hdr.images
    keys   = (i -> i.Key).(images)

    smax = (s, k) -> maximum(x -> getfield(x, k), s)
    slength = (s, k) -> length(unique(x -> getfield(x, k), s))

    vals = [
        (series.Max_No_Slices, :Slice, "slices"),
        (series.Max_No_Echoes, :Echo, "echoes"),
        (series.Max_No_Dynamics, :Dynamic, "dynamics"),
        (series.Max_No_Phases, :Phase, "phases"),
        (series.Max_No_B_Values, :BValue, "bvalues"),
        (series.Max_No_Gradient_Orients, :Grad_Orient, "gradient orientations"),
    ]

    # Check whether data size in series info is consistent with image info
    ws1 = (f, s, i) ->
        "Maximum number of $f not matching: $s (series info) vs. $i (image info)"

    # restrict to mag, real, imag, phas
    ks = filter(k -> k.Type ∈ ("0", "1", "2", "3", "M", "R", "I", "P"), keys)

    if !isempty(ks)
        for (s, k, f) in vals
            i = smax(ks, k)
            i != s && @warn ws1(f, s, i)
        end

        nl = series.No_Label_Types == 0 ? 1 : series.No_Label_Types
        i = slength(ks, :Label_Type)
        i != nl && @warn ws1("label types", series.No_Label_Types, i)
    end

    # Non-exhaustive image info check
    ws2 = (f, m, l) ->
        "Inconsistent number of $f in image info: $m (maximum) vs. $l (length(unique))"

    # image keys
    for (_, k, f) in vals
        m = smax(keys, k)
        l = slength(keys, k)
        l != m && @warn ws2(f, m, l)
    end

    m = smax(keys, :Index)
    l = slength(keys, :Index)
    l != m + 1 && @warn ws2("rec indices", m, l)

    # image info
    slength(images, :Slice_Thickness) != 1 &&
        @warn "Multiple values for `Slice Thickness` found"

    slength(images, :Slice_Gap) != 1 &&
        @warn "Multiple values for `Slice Gap` found"

    slength(images, :Pixel_Spacing) != 1 &&
        @warn "Multiple values for `Pixel Spacing` found"

    slength(images, :Pixel_Size) != 1 &&
        error("Multiple values for `Pixel Size` found")

    slength(images, :Resolution_X) != 1 &&
        error("Multiple values for `Resolution X` found")

    slength(images, :Resolution_Y) != 1 &&
        error("Multiple values for `Resolution Y` found")

    return nothing
end


_parse(type, str) = parse(type, str)

_parse(::Type{String}, str::AbstractString) = String(str)
_parse(::Type{Time}, str::AbstractString) = Time(str, dateformat"H:M:S")
_parse(::Type{Date}, str::AbstractString) =
    try
        Date(str, dateformat"y.m.d")
    catch
        Date(str, dateformat"d.m.y")
    end

_parse(::Type{Bool}, str::AbstractString) =
    str[1] == '0' || uppercase(str[1]) == 'N' ? false :
    str[1] == '1' || uppercase(str[1]) == 'Y' ? true :
    parse(Bool, str)

_parse(::Type{T}, str::AbstractString) where {T<:Tuple} =
    Tuple(_parse.(T.types, split(str)))


function Base.show(io::IO, s::T) where {T<:Union{SeriesInfo, ImageInfo, ImageKey}}
    ks = fieldnames(T)

    println(io, string(T))

    for k in ks
        v = getfield(s, k)
        @printf(io, "    %-27s : ", string(k))
        println(io, v)
    end
end


function Base.show(io::IO, s::RecHeader)
    N = length(s.images)

    println(io, typeof(s))
    println(io)
    print(io, s.series)

    if N > 0
        print(io, "\nImage 1/$(N)\n")
        print(io, s.images[1])
    end

    if N > 1
        println(io)
        print(io, "    .\n")
        print(io, "    .\n")
        print(io, "    .\n")
        println(io)
        print(io, "Image $(N)/$(N)\n")
        print(io, s.images[N])
    end
end


end # module
