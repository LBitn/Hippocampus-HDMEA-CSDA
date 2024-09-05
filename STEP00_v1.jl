# ----------------------------------------------------------------------------------------- #
#=
    Module STEP00_v1
        ] add Dates JLD2 HDF5 InteractiveUtils Suppressor Primes StatsBase Plots Measures
        update
        precompile
        build
=#
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
module STEP00_v1
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Required native packages
# ----------------------------------------------------------------------------------------- #
using Dates
using JLD2
using HDF5
using InteractiveUtils
using Suppressor
using Primes
using StatsBase
using Plots
using Measures
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Functions to be exported
# ----------------------------------------------------------------------------------------- #
export FindDirsFiles
export GetVarsHDF5
export GetChunkSize
export OneSegment
export Digital2Analogue
export SupThr
export UniqueCount
export STDΔV
export PatchEmpties
export Zplot
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Functions
# ----------------------------------------------------------------------------------------- #
"""
    FindDirsFiles( start::String, word::String ) → D::Vector{String}, F::Vector{String}
        Search for directories and files, whose absolute path contains the string "word".
        It begins from the "start" path.
"""
function FindDirsFiles( start::AbstractString, word::AbstractString )
    matching_paths = Vector{ String }( );
    function SearchDir( dir::AbstractString )
        entries = readdir( dir );
        for entry in entries
            entry_path = joinpath( dir, entry );
            if isdir( entry_path ) && !startswith( entry, "." ) # Exclude hidden folders
                if contains( entry, word ) # Case sensitive search for directories
                    push!( matching_paths, entry_path );
                end
                SearchDir( entry_path );
            elseif contains( entry, word ) # Case sensitive search for files
                push!( matching_paths, entry_path );
            end
        end
    end
    SearchDir( start );
    D = matching_paths[ isdir.( matching_paths ) ];
    F = matching_paths[ isfile.( matching_paths ) ];
    return D, F
end

"""
    GetGroupsHDF5( file::HDF5.File, g::AbstractString ) → Groups::Vector{String}
        Extract the groups contained in an HDF5 open file
"""
function GetGroupsHDF5( file::HDF5.File, g::String )
    return [ string( g, "/", key ) for key in keys( file[ g ] ) ];
end

"""
    GetAttrHDF5( file::HDF5.File, g::AbstractString ) → Attributes::Vector{String}
        Extract the attributes contained in an HDF5 open file
"""
function GetAttrHDF5( file::HDF5.File, g::String )
    return [ string( g, "/", key ) for key in keys( attributes( file[ g ] ) ) ];
end

"""
    HDF5Content( file::HDF5.File ) → Groups, Attributes
        Extract Groups and Attributes contained in an open HDF5 file
        # Custom
        using GetGroupsHDF5, GetAttrHDF5
"""
function HDF5Content( file::HDF5.File )
    Groups = keys( file );
    AllGroups = Set( );
    AllAttributes = Set( keys( attributes( file ) ) );
    while !isempty( Groups )
        NewGroups = Set( );
        for g in Groups
            if isa( file[ g ], HDF5.Dataset )
                push!( AllGroups, g );
            else
                union!( NewGroups, GetGroupsHDF5( file, g ) );
                union!( AllAttributes, GetAttrHDF5( file, g ) );
            end
        end
        Groups = collect( NewGroups );
        union!( AllGroups, Groups );
    end
    Groups = sort( unique( collect( AllGroups ) ) );
    Attributes = reverse( sort( unique( collect( AllAttributes ) ) ) );
    return Groups, Attributes
end

"""
    ExtractRawDataset( file::HDF5.File, GroupsHDF5::Vector{ Any } ) → ΔΣ::String, nΔΣ::Vector{String}
        Detects which is the largest dataset contained in the HDF5 file ( ΔΣ::String ).
        It must be marked to avoid loading it or reading it into memory in subsequent steps.
"""
function ExtractRawDataset( file::HDF5.File, GroupsHDF5::Vector{Any} )
    Types = Vector{ String }( undef, length( GroupsHDF5 ) );
    for g = 1:length( GroupsHDF5 )
        Types[ g ] = string( typeof( file[ GroupsHDF5[ g ] ] ) );
    end
    AllDataSets = GroupsHDF5[ Types .== "HDF5.Dataset" ];
    aux = zeros( Int, length( AllDataSets ) );
    for i = 1:length( AllDataSets )
        aux[ i ] = length( file[ AllDataSets[ i ] ] );
    end
    ΔΣ = AllDataSets[ aux .== maximum( aux ) ][ 1 ];
    nΔΣ = AllDataSets[ aux .!= maximum( aux ) ];
    aux = aux[ aux .!= maximum( aux ) ];
    nΔΣ = nΔΣ[ aux .!= 0 ];
    return ΔΣ, nΔΣ
end

"""
    ExtractValues( file::HDF5.File,  AllAttributes::Vector{ String }, nΔΣ::Vector{ Any } ) → D::Dict
        Extracts the values in float or integer that are inside the HDF5, except for the dataset with the largest size.
        Of the latter ( ΔΣ ), only the path is kept for future use.
        # Native
        using HDF5
"""
function ExtractValues( filename::HDF5.File,  AllAttributes::Vector{ String }, nΔΣ::Vector{ Any } )
    D = Dict( );
    e = [ ];
    for g in nΔΣ
        try
            D[ g ] = Float64.( read( filename[ g ] ) );
        catch e
            D[ g ] = read( filename[ g ] );
        end
    end
    for g in keys( D )
        if length( D[ g ] ) .== 1
            D[ g ] = D[ g ][ 1 ];
        end
        try
            D[ g ] = Int64( D[ g ] );
        catch e
        end
    end
    for g in AllAttributes
        aux00 = basename( g );
        aux01 = dirname( g );
        if isempty( aux01 )
            D[ g ] = read_attribute( filename, aux00 );
        else
            D[ g ] = read_attribute( filename[ aux01 ], aux00 );
        end
    end
    return D
end
"""
    Abs2RelPath( D::Dict ) → ND::Dict
        Transforms the key entries from dictionary D from absolute paths to relative paths.
"""
function Abs2RelPath( D::Dict )
    ND = Dict( );
    K = keys( D );
    for k in K
        ND[ basename( k ) ] = D[ k ];
    end
    return ND
end

"""
    ExpDate2Str( Variables::Dict ) → Variables::Dict
        If an "ExperimentDateTime" key exists in the "Variables" dictionary, it is converted to a human-readable format.
        # Native
        using Dates
"""
function ExpDate2Str( Variables::Dict )
    ExperimentDateTime = Variables[ "ExperimentDateTime" ];
    Dt = split( ExperimentDateTime, ":" );
    Dt[ end ] = string( round( Int, parse( Float64, replace( Dt[ end ], r"[A-Z]" => "" ) ) ) );
    newDt = String( "" );
    for i in Dt
        newDt = string( newDt, ":", i );
    end
    newDt = newDt[ 2:end ];
    ExperimentDateTime = Dates.DateTime( newDt );
    Variables[ "ExperimentDateTime" ] = string( ExperimentDateTime );
    return Variables
end

"""
    ExpSett2Dict( Variables::Dict ) → Variables::Dict
        If an "ExperimentSettings" key exists in the "Variables" dictionary, it extracts the values and places them at the same level as the Variables dictionary.
"""
function ExpSett2Dict( Variables::Dict )
    ExperimentSettings = Variables[ "ExperimentSettings" ];
    t = split( ExperimentSettings, "\r\n" );
    t = replace.( t, r"  |{|}|\x22" => "" );
    t = [ i for i in t if !isempty( eachmatch( r"[a-z]", i ) ) ];
    t = split.( t, ": " );
    D = Dict( );
    for i in t
        if i[ 2 ] != ""
            aux = replace( i[ 2 ], r",|\[|\]|\[\]| " => "" );
            if aux != ""
                try
                    aux = parse( Float64, aux );
                catch e
                end
                D[ i[ 1 ] ] = aux;
            end
        end
    end
    delete!( Variables, "ExperimentSettings" );
    Variables = merge( Variables, D );
    return Variables
end

"""
    CleanDictionary( D::Dict ) → D::Dict
        Removes empty entries from a dictionary.
"""
function CleanDictionary( D::Dict )
    delete!( D, String.( keys( D ) )[ values( D ) .== "null" ] );
    delete!( D, [ k for k in keys( D ) if isempty( get( D, k, [ ] ) ) ] );
    return D
end

"""
    right_now( ) → time_stamp::String
        Creates a label with the time (date and time) at which the function is called.
        # Native
        using Dates
"""
function right_now( )
    return Dates.format( round( Dates.DateTime( Dates.now( ) ), Dates.Minute( 15 ) ), Dates.RFC1123Format )
end



"""
    GetVarsHDF5( FILEBRW::String ) → Variables::Dict
        Reads the contents of an HDF5 file. Creates a new directory with the file name.
        Saves the information in a .jld2 file in the same directory.
        Creates a "Variables" variable in the workspace for further manipulation.
        Creates a text report with information from the file.
        # Custom
        using HDF5Content, ExtractRawDataset, ExtractValues, Abs2RelPath, ExpSett2Dict, ExpDate2Str, CleanDictionary
        # Native
        using JLD2, HDF5, InteractiveUtils, Suppressor
"""
function GetVarsHDF5( FILEBRW::String )
    BRW = h5open( FILEBRW, "r" ); # Open the HDF5 file ( not read, DO NOT READ )
    GroupsHDF5, AttrHDF5 = HDF5Content( BRW ); # Extract Groups and Attributes
    ΔΣ, nΔΣ = ExtractRawDataset( BRW, GroupsHDF5 ); # Detect the RAW dataset
    Variables = ExtractValues( BRW,  AttrHDF5, nΔΣ ); # Extract all values except the RAW dataset
    Variables = Abs2RelPath( Variables ); # Keys from Absolute Path to Relative Path
    if ( "ExperimentSettings" in keys( Variables ) )
        Variables = ExpSett2Dict( Variables );
        Variables = ExpDate2Str( Variables );
    end
    Variables = CleanDictionary( Variables );
    dset = BRW[ ΔΣ ];
    filesize = round( ( stat( FILEBRW ).size ) / ( 1024 ^ 3 ), digits = 2 );
    DSETSIZE = deepcopy( filesize );
    try
        DSETSIZE = ( sizeof( dset[ 1 ] ) * size( dset )[ 1 ] ) / ( 1024 ^ 3 );
    catch e
        a, b = size( dset );
        onenumber = h5read( FILEBRW, ΔΣ, ( 1, 1 ) );
        DSETSIZE = ( sizeof( onenumber ) * a * b ) / ( ( 1024 ) ^ 3 );
    end
    Variables[ "RAW" ] = ΔΣ;
    Variables[ "BRWNAME" ] = BRW.filename;
    Variables[ "DSETSIZE" ] = DSETSIZE;
    try
        Variables[ "nChs" ] = length( Variables[ "Chs" ] );
    catch e
        Variables[ "nChs" ] = 4096;
    end
    if !( "NRecFrames" in keys( Variables ) )
        a, b = size( dset );
        if a == 4096
            Variables[ "NRecFrames" ] = b;
        else
            Variables[ "NRecFrames" ] = size( dset, 1 ) / Variables[ "nChs" ];
        end
    end
    BRWTIME = round( ( Variables[ "NRecFrames" ] / Variables[ "SamplingRate" ] ), digits = 3 );
    Variables[ "BRWTIME" ] = BRWTIME;
    PATHBRWs = dirname( FILEBRW );
    PATHMAIN = joinpath( dirname( PATHBRWs ), split( basename( FILEBRW ), "." )[ 1 ] );
    PATHINFO = joinpath( PATHMAIN, "Info" ); mkpath( PATHINFO );
    cd( PATHMAIN )
    FILEVARIABLES = joinpath( PATHINFO, "Variables.jld2" );
    jldsave( FILEVARIABLES; Variables );
    println( "# -------------------------------------------------- Report -------------------------------------------------- #" );
    println( "File: ", replace( BRW.filename, homedir( ) => "~" ) );
    println( "Description: ", Variables[ "Description" ] );
    println( "HDF5 file size: $filesize GB, corresponding to $BRWTIME seconds" );
    println( "Date of Analysis: ", right_now( ), " by ", basename( homedir( ) ) );
    println( "You are now working on the new main path: ", PATHMAIN );
    println( "# ------------------------------------------------------------------------------------------------------------ #" );
    output = @capture_out begin
        println( "# -------------------------------------------------- Report -------------------------------------------------- #" );
        println( "File: ", replace( BRW.filename, homedir( ) => "~" ) );
        println( "Description: ", Variables[ "Description" ] );
        println( "HDF5 file size: $filesize GB, corresponding to $BRWTIME seconds" );
        println( "Date of Analysis: ", right_now( ), " by ", basename( homedir( ) ) );
        println( "With:\n", string( basename( homedir( ) ), "@", gethostname( ) ) );
        versioninfo( )
        println( "# ------------------------------------------------------------------------------------------------------------ #" );
    end
    FILEINFO = joinpath( PATHINFO, "Reporte.txt" );
    write( FILEINFO, output )
    close( BRW )
    return Variables
end

"""
    GetChunkSize( Variables::Dict, MaxGB::Real = 0.5, m::Int = 3, M::Int = 500 )
        Calculates the number of segments in which to cut the dataset, whose disk space does not exceed "MaxGB".
        If it is not possible to determine a suitable number of segments, resort to the default value, 4 seconds per segment.
        It also discards the results with more than 500 segments or less than 3 segments.
        These limits are optional and can be modified in the "m" and "M" input.
        # Native
        using Primes
"""
function GetChunkSize( Variables::Dict, MaxGB::Real = 0.5, m::Int = 3, M::Int = 500 )
    NRecFrames = Variables[ "NRecFrames" ];
    DSETSIZE = Variables[ "DSETSIZE" ];
    nChs = Variables[ "nChs" ];
    SamplingRate = Variables[ "SamplingRate" ];
    onenumber = DSETSIZE / ( NRecFrames * nChs );
    while isprime( NRecFrames );
        NRecFrames = NRecFrames - 1;
        DSETSIZE = DSETSIZE - onenumber;
    end
    divs = divisors( NRecFrames );
    divs = divs[ divs .>= m .&& divs .<= M ];
    aux = DSETSIZE ./ divs;
    aux = aux[ aux .<= MaxGB ];
    finalsize = 0;
    σ = 0;
    if isempty( aux )
        println( "No optimal segment size has been found for the “MaxGB” value. We will use the default segment length of 4 seconds." );
        σ = floor( Int, NRecFrames / SamplingRate * 4 ); # 4 secs
        nfrs = floor( Int, NRecFrames / σ );
        NRecFrames = Int( nfrs * σ );
        finalsize = nfrs * nChs * onenumber;
    else
        finalsize = maximum( aux );
        σ = Int( DSETSIZE  / finalsize );
    end
    finaltime = ( ( NRecFrames / σ ) / SamplingRate );  # sec
    fs = round( finalsize, digits = 3 );
    ft = round( finaltime, digits = 3 );
    println( "$σ segments of $fs GB and $ft seconds each" );
    return σ
end

"""
    OneSegment( Variables::Dict, n::Int, N::Int ) -> BIN::Array{ UInt16 }
        Cut the n segment ( n-th from N ) from the dataset.
        The result array has n-channels X n-frames form.
        # Native
        using HDF5
"""
function OneSegment( Variables::Dict, n::Int, N::Int )
    RAW = h5open( Variables[ "BRWNAME" ], "r" )[ Variables[ "RAW" ] ];
    NRecFrames = Variables[ "NRecFrames" ];
    nfrs = floor( Int, ( NRecFrames / N ) );
    BIN = Array{ UInt16 }( undef, Variables[ "nChs" ], nfrs );
    nChs = Variables[ "nChs" ];
    nchs = [ ];
    nFrs = [ ];
    try
        nchs, nFrs = size( RAW );
    catch e
        nchs = 1;
        nFrs = length( RAW );
    end
    if nchs == nChs # the dataset has an n x m form. Older brw versions
        fr0 = ( ( ( n - 1 ) * nfrs ) + 1 );
        frN = n*nfrs;
        BIN = RAW[ ( 1 : nChs ), ( fr0 : frN ) ];
    elseif nFrs == NRecFrames*nChs # the dataset has a vector form ( 1 x m ). Newest brw versions
        for frame in ( ( n - 1 ) * nfrs + 1 ) : ( nfrs * n )
            f = frame - ( nfrs * ( n - 1 ) );
            BIN[ :, f ] = RAW[ ( ( frame - 1 ) * nChs + 1 ) : ( nChs * frame ) ];
        end
    end
    close( RAW )
    return BIN
end

"""
    Digital2Analogue( Variables::Dict, DigitalValue::Matrix{UInt16} ) → BIN::Matrix{Float64}
        Conversion of ΔΣ data extracted from the brw file to voltage values (μV) for Matrix format acording to the equation:
        Voltage = ( DigitalValue + ADCCountsToMV ) * MVOffset
        DigitalValue == ΔΣ data in Matrix{ UInt16 } form
"""
function Digital2Analogue( Variables::Dict, DigitalValue::Matrix{ UInt16 } )
    SignalInversion = Variables[ "SignalInversion" ];
    MinVolt = Variables[ "MinVolt" ];
    MaxVolt = Variables[ "MaxVolt" ];
    BitDepth = Variables[ "BitDepth" ];
    MVOffset = SignalInversion * MinVolt;
    ADCCountsToMV = ( SignalInversion * ( MaxVolt - MinVolt ) ) / ( 2 ^ BitDepth );
    BIN = @. MVOffset + ( DigitalValue * ADCCountsToMV );
    return BIN
end

"""
    SupThr( Data::Array, Thr::Real ) → Cols::Vector, Rows::Vector
        Find values above Thr and below -Thr. Organized on Columns and Rows
"""
function SupThr( Data::Array, Thr::Real )
 ST = findall( abs.( Data ) .>= Thr );
    AllCols = getindex.( ST, [ 1 ] );
    AllRows = getindex.( ST, [ 2 ] );
    Rows = [ ];
    Cols = [ ];
    for col in sort( unique( AllCols ) )
        push!( Rows, AllRows[ AllCols .== col ] );
        push!( Cols, col )
    end
    return Cols, Rows
end

"""
    UniqueCount( data::Array ) → Count::Vec{Int64}
        Counts the number of unique values for each row of an array
        Sort of Cardinality measure
"""
function UniqueCount( data::Array )
    X, Y = size( data );
    Count = Array{ Int64 }( undef, X );
    [ Count[ x ] = length( unique( round.( data[ x, : ], digits = 2 ) ) ) for x in 1:X ];
    return Count
end

"""
    STDΔV( Variables::Dict, BIN::VecOrMat, ΔT::Real ) → STD::Vector{Float64}
        ΔT in msecs
        # Custom
        using ms2frs
        # Native
        using StatsBase
"""
function STDΔV( Variables::Dict, BIN::VecOrMat, ΔT::Real = 250 )
    BIN = Float64.( BIN );
    ΔT = ms2frs( ΔT, Variables );
    STD = vec( std( ( circshift( BIN, ( 0, ΔT ) ) - BIN ), dims = 2 ) );
    return STD
end

"""
    ms2frs( time::Real, SamplingRate::Real OR Variables::Dict ) → frs::Int
        Converts msecs to frs in accordance with the sampling rate of the file
"""
function ms2frs( time::Real, SamplingRate::Real )
    return ceil( Int, ( time * SamplingRate ) / 1000 );
end
function ms2frs( time::Real, Variables::Dict )
    SamplingRate = Variables[ "SamplingRate" ];
    return ceil( Int, ( time * SamplingRate ) / 1000 );
end

"""
    PatchEmpties( aux::Vector, Empties::Vector = [ ] ) → aux::Vector
        Replaces the aux vector values in the Empties vector positions with random non-Empties aux values.
        For graphing purposes only.
        # Native
        using StatsBase
"""
function PatchEmpties( aux::Vector, Empties::Vector = [ ] )
    nChs = length( aux );
    NotEmpties = setdiff( 1:nChs, Empties );
    nv = sample( NotEmpties, length( Empties ) );
    aux[ Empties ] = aux[ nv ];
    return aux
end

"""
    Zplot( W::VecOrMat, cm_::Symbol, t::String, c::Real = 0 ) → P::Plot
        Plot a square heatmap, gr backend, cm_ colormap with title "t", cbarlims = c
        If the color bar limits are not defined (c), it will be set to ± ( media + 2*std ) of W
        # Native
        using Plots, Measures, StatsBase
"""
function Zplot( W::VecOrMat, cm_::Symbol = :vik, t::String = "", c::Real = 0 )
    W = vec( W );
    nc = Int( sqrt( length( W ) ) );
    Z = reverse( reshape( W, nc, nc )', dims = 1 );
    if c == 0
        c = median( W ) + ( 2 * std( W ) );
    end
    Plots.gr( );
    P = plot( );
    P = heatmap(
        Z,
        wsize = ( 400, 400 ),
        axis = ( [ ], false ),
        title = t,
        titlefontsize = 10,
        cbarfontsize = 10,
        dpi = 300,
        fillalpha = 1.0,
        right_margins = 5mm,
        clims = ( -c, c ),
        colormap = cm_,
        aspect_ratio = :equal,
        lims = ( 0, nc + 1 ),
        );
    return P
end
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
end # module STEP00_v1
# ----------------------------------------------------------------------------------------- #
