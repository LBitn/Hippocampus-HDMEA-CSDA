# ----------------------------------------------------------------------------------------- #
#=
    Module STEP01_v1
        ] add JLD2 StatsBase Plots Measures
        update
        precompile
        build
=#
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
module STEP01_v1
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Required native packages
# ----------------------------------------------------------------------------------------- #
using JLD2
using StatsBase
using Plots
using Measures
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Functions to be exported
# ----------------------------------------------------------------------------------------- #
export FindDirsFiles
export SearchDir
export LoadDict
export BarPlot
export SupThr
export ReduceArrayDistance
export Neighbours
export ReconstructChannels
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
    SearchDir( path::String, key::String ) → list::Vector{ String }
"""
SearchDir( path::String, key::String ) = filter( x -> endswith( x, key ), readdir( path; join = true ) );

"""
    LoadDict( filename::String ) → D::Dict
        # Native
        using JLD2
"""
function LoadDict( filename::String )
    D = load( filename );
    K = keys( D );
    if length( K ) == 1
        k = collect( K )[ 1 ];
        D = D[ k ];
    end
    return D
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
    BarPlot( W::VecOrMat, fc::Symbol = :royalblue3, t::String = "", xl::String = "", yl::String = "" ) → Plot
        # Native
        using Plots, Measures
"""
function BarPlot( W::VecOrMat, fc::Symbol = :royalblue3, t::String = "", xl::String = "", yl::String = "" )
    P = plot( );
    P = bar( W,
        leg = :none,
        fillcolor = fc,
        lc = :white,
        grid = :none,
        wsize = ( 600, 400 ),
        title = t,
        top_margin = 10mm,
        dpi = 300,
        xlabel = xl,
        ylabel = yl,
        left_margin = 5mm,
        bottom_margin = 5mm,
        );
    return P
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
    ReduceArrayDistance( W::Vector, distance::Int64 ) → G::Array
        Groups contiguous numbers with maximum distance between them “distance”.
"""
function ReduceArrayDistance( W::Vector, distance::Int64 )
    g = [ ];
    W = sort( unique( W ) );
    for w0 in W
        r = collect( w0 :( w0 + distance ) );
        push!( g, W[ W .∈ [ r ] ] );
    end
    G = [ ];
    push!( G, g[ 1 ] );
    for i = 2:length( g )
        if isempty( intersect( g[ i ], G[ end ] ) )
            push!( G, g[ i ] );
        else
            G[ end ] = sort( union( g[ i ], G[ end ] ) );
        end
    end
    return G
end

"""
    Neighbours( C::Int64, d::Int64 ) → A::Array{ Int64 }, v::Vector{ Int64 }
        A = Array( ( d*2 ) + 1, ( d * 2 ) + 1 ),
        v = vec( 2*( ( d * 2 ) + 1 ) - 1 );
        The d-neighborhood is calculated from the channel ( C ) as a center
        A = array where C is the center and is in chip order
        v = same neighboring channels as A but in vector form and without C ( 8 channels for d = 1 )
"""
function Neighbours( center::Int64, d::Int64 )
    Layout = reverse( reshape( collect( 1:4096 ), 64, 64 )', dims = 1 );
    x = findall( Layout .== center )[ ][ 2 ];
    y = findall( Layout .== center )[ ][ 1 ];
    aux = [ ( x - d ),( x + d ), ( y - d ), ( y + d ) ];
    aux[ aux .< 1 ] .= 1;
    aux[ aux .> 64 ] .= 64;
    A = Layout[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    v = vec( A )[ vec( A ) .!= center ];
    return A, sort( v )
end

"""
    Fences( data::Vector ) → LowerFence::Real, HigherFence::Real
        Simple test for Outliers
        # Native
        using StatsBase
"""
function Fences( data::Vector )
    Q1 = quantile( data, 0.25 );
    Q3 = quantile( data, 0.75 );
    IQR = Q3 - Q1;
    LF = Q1 - 1.5*IQR;
    HF = Q3 + 1.5*IQR;
    return LF, HF
end

"""
    ReconstructChannels( data::Array, lim = 10 ) → fictional_channel::Vector
        Creates an vector form an array in the following order: "mean", "median", "n random samples"
        testing the std form the array vs the resultant vector. If the std form the vector is an
        outlier acordingly on the fences test, then the next vector is determined. If is necesary to
        take n random samples, the limit number of times is set by the user
        # Custom
        using Fences
        # Native
        using StatsBase
"""
function ReconstructChannels( data::Array, lim = 10 )
    _, nFrs = size( data );
    STDS = vec( std( data, dims = 2 ) );
    LF, HF = Fences( STDS )
    C = 0;
    fictional_channel = vec( mean( data, dims = 1 ) );
    s = std( fictional_channel );
    test = ( s > LF && s < HF );
    if test == false
        fictional_channel = vec( median( data, dims = 1 ) );
        s = std( fictional_channel );
        test = ( s > LF && s < HF );
        if test == false
            while test == false && C < lim
                fictional_channel = sample( data, nFrs );
                s = std( fictional_channel );
                test = ( s > LF && s < HF );
                C = C + 1;
            end
        end
    end
    return fictional_channel
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
function STDΔV( Variables::Dict{Any, Any}, BIN::VecOrMat, ΔT::Real = 250 )
    BIN = Float64.( BIN );
    ΔT = ms2frs( ΔT, Variables );
    STD = vec( std( ( circshift( BIN, ( 0, ΔT ) ) - BIN ), dims = 2 ) );
    return STD
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
end # module STEP01_v1
# ----------------------------------------------------------------------------------------- #

