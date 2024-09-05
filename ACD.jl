# ----------------------------------------------------------------------------------------- #
#=
    Module ACD
        ] add JLD2 StatsBase BinningAnalysis HistogramThresholding Suppressor Plots Measures
        update
        precompile
        build
=#
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
module ACD
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Required native packages
# ----------------------------------------------------------------------------------------- #
using JLD2
using StatsBase
using BinningAnalysis
using HistogramThresholding
using Suppressor
using Plots
using Measures
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Functions to be exported
# ----------------------------------------------------------------------------------------- #export FindDirsFiles
export LoadDict
export PatchEmpties
export GetSupThr
export GetGroups
export Neighbours
export GetPatches
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
    Thresholding( aux::VecOrMat ) -> t::vec
        # Native
        using BinningAnalysis, StatsBase, HistogramThresholding, Suppressor
        # Custom
        using Donoho, SigmaData, Fences
"""
function Thresholding( aux::VecOrMat )
    W = round.( vec( Float64.( aux ) ), digits = 3 );
    t = zeros( 19 );
    @suppress begin
        xmean, xerror = jackknife( identity, W );
        t[ 1 ] = Donoho( W ) * SigmaData( W );
        t[ 2 ] = mean( W ) + ( 2 * std( W ) );
        t[ 3 ] = mean( W ) - ( 2 * std( W ) );
        t[ 4 ], t[ 5 ] = Fences( W );
        try t[ 6 ] = find_threshold( W, Balanced( ) ); catch e; end
        try t[ 7 ] = find_threshold( W, Entropy( ) ); catch e; end
        try t[ 8 ] = find_threshold( W, Intermodes( ) ); catch e; end
        try t[ 9 ] = find_threshold( W, MinimumError( ) ); catch e; end
        try t[ 10 ] = find_threshold( W, MinimumIntermodes( ) ); catch e; end
        try t[ 11 ] = find_threshold( W, Moments( ) ); catch e; end
        try t[ 12 ] = find_threshold( W, Otsu( ) ); catch e; end
        try t[ 13 ] = find_threshold( W, UnimodalRosin( ) ); catch e; end
        try t[ 14 ] = find_threshold( W, Yen( ) ); catch e; end
        t[ 15 ] = xmean + xerror;
        t[ 16 ] = xmean - xerror;
        t[ 17 ] = percentile( W, 0.33 );
        t[ 18 ] = median( W ) + xerror;
        t[ 19 ] = median( W ) - xerror;
    end
    l, h = extrema( W );
    t = t[ h .>= t .>= l ];
    t = sort( unique( t ) );
    return t
end

"""
    Donoho( x ) = ( median( abs.( x ) ) / 0.6745 )
        # Native
        using StatsBase
"""
Donoho( x ) = ( median( abs.( x ) ) / 0.6745 );

"""
    SigmaData( data ) = sqrt( 2 * log( length( data ) ) ) -> σ::Real
"""
SigmaData( data ) = sqrt( 2 * log( length( data ) ) );

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

"""
    Z0( X::VecOrMat, nChs::Int = 4096 ) -> P::Plot
        # Native
        using Plots
"""
function Z0( X::VecOrMat, nChs::Int = 4096 )
    X = Int.( vec( X ) );
    Z = zeros( Int, nChs );
    n = Int( sqrt( nChs ) );
    Z[ X ] .= Z[ X ] .+ 1;
    Z = Int.( reverse( reshape( Z, n, n )', dims = 1 ) );
    Plots.gr( );
    P = plot( );
    P = heatmap( Z,
        aspect_ratio = :equal,
        wsize = ( 400, 400 ),
        axis = ( [ ], false ),
        dpi = 300,
        margins = 5mm,
        cbar = :none,
        colormap = :greys,
        lims = ( 0, n + 1 )
        );
    return P
end

"""
    GetInfThr( W, upper_limit ) -> r, P
        # Custom
        using Zplot, Thresholding, Z0
        # Native
        using Plots, Measures
"""
function GetInfThr( W, upper_limit )
    P0 = Zplot( W );
    t = Thresholding( W );
    M = repeat( W, 1, length( t ) );
    aux0 = vec( sum( M .<= t', dims = 2 ) );
    P1 = Zplot( aux0, :greys );
    t0 = length( t );
    r = findall( aux0 .>= t0 );
    while length( r ) <= upper_limit && t0 > 2
        t0 = t0 - 1;
        r = findall( aux0 .>= t0 );
    end
    nchs = length( r );
    P2 = Z0( r );
    P = plot(
        P0, P1, P2,
        layout = ( 1, 3 ),
        wsize = ( 800, 400 ),
        title = [""  "\n" ^ 2 *"$nchs Detected" ""],
        titlefont = Plots.font( pointsize = 10, family = "sans-serif" ),
        margin = 5mm,
        cbar = :none
    );
    return r, P
end

"""
    GetSupThr( W, lower_limit::Int = 1500 ) -> r, P
        # Custom
        using Zplot, Thresholding, Z0
        # Native
        using Plots, Measures
"""
function GetSupThr( W, lower_limit::Int = 1500 )
    P0 = Zplot( W );
    t = Thresholding( W );
    M = repeat( W, 1, length( t ) );
    aux0 = vec( sum( M .>= t', dims = 2 ) );
    P1 = Zplot( aux0, :greys );
    t0 = length( t );
    r = findall( aux0 .>= t0 );
    while length( r ) <= lower_limit && t0 > 2
        t0 = t0 - 1;
        r = findall( aux0 .>= t0 );
    end
    nchs = length( r );
    P2 = Z0( r );
    P = plot(
        P0, P1, P2,
        layout = ( 1, 3 ),
        wsize = ( 800, 400 ),
        title = [""  "\n" ^ 2 *"$nchs Detected" ""],
        titlefont = Plots.font( pointsize = 10, family = "sans-serif" ),
        margin = 5mm,
        cbar = :none
    );
    return r, P
end

"""
    GetGroups( W ) -> singles, Gs
        # Custom
        using Neighbours, LoopGroups
"""
function GetGroups( W )
    Gs = [ ];
    for ch in W
        neigs, _ = Neighbours( ch, 1 );
        aux = neigs[ neigs .∈ [ W ] ];
        push!( Gs, aux )
    end
    singles = Gs[ length.( Gs ) .== 1 ];
    singles = vcat( singles... );
    Gs = Gs[ length.( Gs ) .!= 1 ];
    Gs = sort( unique( sort.( Gs ) ) );
    a = length( Gs );
    Gs = LoopGroups( Gs );
    b = length( Gs );
    while a != b
        a = length( Gs );
        Gs = LoopGroups( Gs );
        b = length( Gs );
    end
    return singles, Gs
end

"""
    LoopGroups( Gs ) -> GS
"""
function LoopGroups( Gs )
    GS = [ ];
    i = 1
    while i <= length( Gs )
        NewGroup = [ ];
        gs = Gs[ i ];
        push!( NewGroup, gs );
        for j = 2:length( Gs )
            push!( NewGroup, intersect( gs, Gs[ j ] ) );
        end
        aux = length.( NewGroup ) .!= 0;
        if length( aux ) > 1
            g = unique( vcat( Gs[ aux ]... ) );
        else
            g = Gs[ aux ];
        end
        push!( GS, vcat( g... ) );
        deleteat!( Gs, aux );
    end
    return GS
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
    GetPatches( G::Vector{Any}, thr::Int = 4 ) -> chs, P
        # Custom
        using Z0
        # Native
        using Plots
"""
function GetPatches( G::Vector{Any}, thr::Int = 4 )
    g = G[ length.( G ) .>= thr ];
    chs = sort( vcat( g... ) );
    nchs = length( chs );
    Plots.gr( );
    P = plot( );
    P = Z0( chs );
    P = plot!( P,
        title = "\n" ^ 2 *"$nchs Selected, thr = $thr",
        titlefont = Plots.font( pointsize = 10, family = "sans-serif" )
    );
    return chs, P
end
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
end # module ACD
# ----------------------------------------------------------------------------------------- #



















