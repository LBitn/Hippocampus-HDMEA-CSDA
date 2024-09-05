# ----------------------------------------------------------------------------------------- #
#=
    Module STEP02_v1
        ] add JLD2
        update
        precompile
        build
=#
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
module STEP02_v1
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Required native packages
# ----------------------------------------------------------------------------------------- #
using JLD2
using StatsBase
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Functions to be exported
# ----------------------------------------------------------------------------------------- #
export FindDirsFiles
export LoadDict
export SearchDir
export CSDA
export GetCentersOfMass
export Trajectories
export FixingGaps
export StartStop
export DonohoMatrix3D
export SigmaData
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Constants
# ----------------------------------------------------------------------------------------- #
GaussianKernel = [
    0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
    0.00002292	0.00078634	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
    0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
    0.00038771	0.01330373	0.11098164	0.22508352	0.11098164	0.01330373	0.00038771
    0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
    0.00002292	0.00078633	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
    0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
    ];
# The Laplace-Lindenberg operator
Laplacian01 = [ [ 0 1 0 ]; [ 1 -4 1 ]; [ 0 1 0 ] ];
Laplacian02 = [ [ ( 1 / 2 ) 0 ( 1 / 2 ) ]; [ 0 -2 0 ]; [ ( 1 / 2 ) 0 ( 1 / 2 ) ] ];
LaplacianKernel = ( ( 1 - ( 1 / 3 ) ) * Laplacian01 ) + ( ( 1 / 3 ) * Laplacian02 );
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
    SearchDir( path::String, key::String ) → list::Vector{ String }
"""
SearchDir( path::String, key::String ) = filter( x -> endswith( x, key ), readdir( path; join = true ) );


"""
    UnNormGauss( x, σ = 3 ) → exp( ( -x * x ) / ( 2 * σ ) )
"""
function UnNormGauss( x, σ::Real = 3 )
    return exp( ( -x * x ) / ( 2 * σ ) );
end

"""
    GaussSmoothTemporal( Data::Array, σ = 3 ) → GS::Array
        A temporal Gaussian smoothing.
        This is essentially a low-pass filter.
        It depends implicitly on the sampling frequency.
        σ is measured in pixels, it is the standard deviation of our kernel.
        The α of our window will be 3 * σ.
"""
function GaussSmoothTemporal( Data::Array, σ::Real = 3 )
    α = ceil( Int, σ * 3 );
    β = ones( α );
    GS = zeros( size( Data ) );
    aux = vcat( β * Data[ 1 ], Data, β * Data[ end ] );
    kernel = map( x -> UnNormGauss( x, σ ), collect( -α:α ) );
    kernel = kernel / ( sum( kernel ) );
    for t = ( α + 1 ):( length( Data ) + α )
        GS[ t - α ] = sum( aux[ ( t - α ):( t + α ) ] .* kernel );
    end
    return GS
end

"""
    DiscreteLaplacian( Data::Array ) → DL::Array
"""
function DiscreteLaplacian( Data::Array )
    ( J, K ) = size( Data );
    I = reshape( Data[ 1, : ], ( 1, K ) );
    R = reshape( Data[ end, : ], ( 1, K ) );
    # Padding with copy of data
    Data = vcat( I, Data, R );
    Data = hcat( Data[ :, 1 ], Data, Data[ :, end ] );
    J, K = size( Data );
    Fx = Array{ Float64 }( undef, 3, 3 );
    DL = zeros( size( Data ) );
    #= calculate the CSD by applying the Laplacian Kernel to each cell plus its 8-Neighborhood
    and sum all the results as th new cell value =#
    for j = ( 2:( J - 1 ) ), k = ( 2:( K - 1 ) )
        Fx = Data[ ( j - 1 ):( j + 1 ), ( k - 1 ):( k + 1 ) ];
        DL[ j, k ] = sum( LaplacianKernel .* Fx );
    end
    # Crop the borders
    DL = DL[ 2:( end - 1 ), 2:( end - 1 ) ];
    return DL
end

"""
    GaussianBlur( Data::Array ) → GB::Array
"""
function GaussianBlur( Data::Array )
    GB = zeros( size( Data ) );
    ( J, K ) = size( Data );
    # Padding with copy of data
    A = reshape( Data[ 1, : ], ( 1, K ) );
    B = reshape( Data[ end, : ], ( 1, K ) );
    C = vcat( A, A, A );
    D = vcat( B, B, B );
    Data = vcat( C, Data, D );
    for j = 1:3
        Data = hcat( Data[ :, 1 ], Data, Data[ :, end ] );
    end
    for j = ( 4:( J + 3 ) ), k = ( 4:( K + 3 ) )
        bin = Data[ ( j - 3 ):( j + 3 ), ( k - 3 ):( k + 3 ) ];
        GB[ ( j - 3 ), ( k - 3 ) ] = sum( GaussianKernel .* bin );
    end
    return GB
end

"""
    CSDA( Data::Array, σ = 3 ) → ∇::Array
        # Custom
        using GaussSmoothTemporal, GaussianBlur, DiscreteLaplacian
"""
function CSDA( Data::Array, σ::Real = 3 )
    nChs, nFrs = size( Data );
    side = Int( sqrt( nChs ) );
    Data = reshape( Data, side, side, nFrs );
    ( J, K, L )  = size( Data );
    # We apply a Temporal Gaussian smoothing ( this greatly affects the animations )
    BINSMOOTH = zeros( J, K, L );
    for j = 1:J, k = 1:K
        BINSMOOTH[ j, k, : ] = GaussSmoothTemporal( vec( Data[ j, k, : ] ), σ );
    end
    ∇ = zeros( J, K, L );
    # We spatially smooth the raw data with a two-dimensional Gaussian filter.
    # Later we obtain the dCSD.
    for l = 1:L
        ∇[ :, :, l ] = DiscreteLaplacian( GaussianBlur( BINSMOOTH[ :, :, l ] ) );
    end
    ∇ = -1 .* ∇;
    return ∇
end

"""
    Donoho( x ) = ( median( abs.( x ) ) / 0.6745 ) → |R|::Real
        using StatsBase
"""
Donoho( x ) = abs( median( abs.( x ) ) / 0.6745 );

"""
    DonohoMatrix3D( Data::Array{ Float64, 3 } ) → R::Matrix{ Float64 }
        using Donoho
"""
function DonohoMatrix3D( Data::Array{ Float64, 3 } )
    a, b, c = size( Data );
    R = Array{ Any }( undef, a, b ); fill!( R, [ ] );
    [ R[ k, j ] = Donoho( Data[ k, j, : ] ) for k in 1:a, j in 1:b ];
    return R
end

"""
    SigmaData( data ) = sqrt( 2 * log( length( data ) ) ) → R::Real
"""
SigmaData( data ) = sqrt( 2 * log( length( data ) ) );

"""
    GetCentersOfMass( CSD::Array{ Float64, 3 }, minchannels::Int, ϵ::Union{ Array, Real } = 0 ) → CMN::Dict{Int64, Array}, CMP::Dict{Int64, Array}
        using MassCenters
        using StatsBase
"""
function GetCentersOfMass( CSD::Array{ Float64, 3 }, minchannels::Int, ϵ::Union{ Array, Real
} = 0 )
    ( rows, cols, tN ) = size( CSD );
    CMP = Dict{ Int, Array }( ); # Positive Mass Centers
    CMN = Dict{ Int, Array }( ); # Negative Mass Centers
    if ϵ == 0
        ϵ = abs.( std( CSD, dims = 3 )[ :, :, 1 ] ) .* 3; # obtain the thr if not given
    end
    for t in 1:tN
        NegChs = Array{ Int16 }[ ];
        PosChs = Array{ Int16 }[ ];
        for row in 1:rows, col in 1:cols
            thr = [ ]; # threshold determination
            if typeof( ϵ ) <: VecOrMat
                thr = abs( ϵ[ row, col ] );
                elseif typeof( ϵ ) <: Real
                thr = ϵ;
            end
            if ( CSD[ row, col, t ] <= ( -1 * thr ) )
                push!( NegChs, [ row, col ] ); # Negative Supra Threshold
                elseif ( CSD[ row, col, t ] >= thr )
                push!( PosChs, [ row, col ] ); # Positive Supra Threshold
            end
        end
        CMN[ t ] = MassCenters( CSD, t, NegChs, minchannels ); # declare sinks
        CMP[ t ] = MassCenters( CSD, t, PosChs, minchannels ); # declare sources
    end
    return CMN, CMP
end

"""
    MassCenters( CSD::Array{ Float64, 3 }, t::Int, channels::Vector, minchannels::Int ) → centers::Matrix{Int64}
        centers = [ x y Ω ]; Cartesian coordinates and Weight
        using DisjointComponents
"""
function MassCenters( CSD::Array{Float64, 3}, t::Int, channels::Vector, minchannels::Int )
    DCs = DisjointComponents( channels );
    centers = [ 0 0 0 ];
    for component in DCs
        μ = length( component );
        if μ >= minchannels
            Ω = 0.00;
            x = 0.00;
            y = 0.00;
            for channel in component
                ω = CSD[ channel[ 1 ], channel[ 2 ], t ];
                Ω += ω;
                x += channel[ 2 ] * ω;
                y += channel[ 1 ] * ω;
            end
            x /= Ω;
            y /= Ω;
            centers = vcat( centers, [ x y Ω ] );
        end
    end
    centers = centers[ ( 2 : end ), : ]; # remove first [ 0 0 0 ] entry
    return centers
end

"""
    DisjointComponents( CartesianChannels::Vector{ Array } ) → componentes::Set{ Any }
        using EightNeigh
"""
function DisjointComponents( CartesianChannels::Vector{ Array{ Int16 } } )
    temp = copy( CartesianChannels );
    components = Set{ Any }( );
    while ( length( temp ) != 0 )
        aux = pop!( temp ); # Removes the LAST item from the list
        aux1 = Array{ Int64 }[ ];
        aux2 = Array{ Int64 }[ ];
        push!( aux1, aux ); # Puts items at the END of the list
        push!( aux2, aux );
        depth = 0;
        while ( length( aux1 ) != 0 )
            one = pop!( aux1 );
            for neigh in EightNeigh( one )
                if in( neigh, temp )
                    deleteat!( temp, indexin( Any[ neigh ], temp ) );
                    push!( aux1, neigh );
                    push!( aux2, neigh );
                end
            end
        end
        push!( components, aux2 );
    end
    return components
end

"""
    EightNeigh( cartesian_channel::Vector{ Int } ) → neighbors::Set{ Any }
        The eight-neighborhood of a point on a square grid.
"""
function EightNeigh( cartesian_channel::Vector{ Int } )
    j, k = cartesian_channel;
    result = Set{ Array{ Int64, 1 } }( );
    push!( result, [ j - 1, k - 1 ] );
    push!( result, [ j - 1, k ] );
    push!( result, [ j - 1, k + 1 ] );
    push!( result, [ j, k - 1 ] );
    push!( result, [ j, k + 1 ] );
    push!( result, [ j + 1, k - 1 ] );
    push!( result, [ j + 1, k ] );
    push!( result, [ j + 1, k + 1 ] );
    neighbors = Set{ Array{ Int64, 1 } }( );
    for p in result
        if !( p[ 1 ] == 0 || p[ 2 ] == 0 || p[ 1 ] == 65 || p[ 2 ] == 65 )
            push!( neighbors, p );
        end
    end
    return neighbors
end

"""
    DistanceVectors( v, M ) → distance
"""
function DistanceVectors( v, M )
    distance = sqrt.( ( ( v[ 1 ] .- M[ :, 1 ] ) .^ 2 ) .+ ( ( v[ 2 ] .- M[ :, 2 ] ) .^ 2 ) );
    return distance
end

"""
    DistanceCoords( x, y ) → result
"""
function DistanceCoords( x, y )
    result = sqrt( ( ( x[ 1 ] - y[ 1 ] ) ^ 2 ) + ( ( x[ 2 ] - y[ 2 ] ) ^ 2 ) );
    return result
end

"""
    Trajectories( CM, tol_dist, tol_time, Tmin, min_weight )
        using WeightSelection
        using AuxiliaryTrajectories
"""
function Trajectories( CM, tol_dist, tol_time, Tmin, min_weight )
    cms = copy( CM );
    tFrs = length( cms );
    for k in 1:5
        cms[ tFrs + k ] = Array{ Float64 }( undef, 0, 3 );
    end
    UltimoFr = false;
    allTjs = Dict{ Integer, Array{ Any } }( );
    nTs = 1;
    for tiempo in 1:tFrs
        if tiempo == tFrs
            UltimoFr = true;
        end
        ConjEnFr = WeightSelection( cms[ tiempo ], min_weight );
        NumConj = size( ConjEnFr )[ 1 ];
        if NumConj > 0
            for j in 1:NumConj
                if tiempo < tFrs
                    UltimoFr = false;
                end
                tprim = tiempo;
                More = true;
                ConjConFr = [ transpose( ConjEnFr[ j, : ] ) tiempo ];
                Tj = ConjConFr;
                t_aux = 0;
                while More == true && t_aux <= tol_time
                    if tprim >= tFrs
                        UltimoFr = true;
                    end
                    if tprim < tFrs
                        ConjEnSigFr = WeightSelection( cms[ tprim + 1 ], min_weight );
                        NumConjSig = size( ConjEnSigFr )[ 1 ];
                    else
                        NumConjSig = 0;
                        t_aux = tol_time;
                    end
                    if NumConjSig > 0
                        ( Tj, cms, More, t_aux, allTjs, nTs ) =
                            AuxiliaryTrajectories(
                                Tj, cms, tprim, tol_dist, More, t_aux, tol_time, allTjs,
                                nTs, Tmin );
                            if UltimoFr == true
                                if size( Tj )[ 1 ] > Tmin
                                    allTjs[ nTs ] = Tj;
                                    nTs += 1;
                                    More = false;
                                end
                            end
                            tprim += 1;
                    else
                        if t_aux == tol_time
                            if size( Tj )[ 1 ] > Tmin
                                allTjs[ nTs ] = Tj;
                                nTs += 1;
                                More = false;
                            else
                                More = false;
                            end
                        end
                        t_aux += 1;
                        tprim += 1;
                    end
                end
            end
        end
    end
    return allTjs
end

"""
    WeightSelection( Fr, min_weight ) -> good_ones::Int
"""
function WeightSelection( Fr, min_weight )
    NF = size( Fr )[ 1 ];
    if NF > 0
        Ws = Fr[ :, 3 ];
        GoodOnes = findall( abs.( Ws ) .> min_weight );
        good_ones = Fr[ GoodOnes, : ];
    else
        good_ones = Fr;
    end
    return good_ones
end

"""
    AuxiliaryTrajectories( Tj, cms, tiempo, tol_dist, More, t_aux, tol_time, allTjs, nTs, Tmin ) -> Chain, cms, More, t_aux, allTjs, nTs
        using Velocities, DistanceVectors
"""
function AuxiliaryTrajectories( Tj, cms, tiempo, tol_dist, More, t_aux, tol_time, allTjs, nTs,
Tmin )
    if size( Tj )[ 1 ] < 5
        tol_dist = tol_dist;
    else
        Locs = Tj[ ( end - 4 ) : ( end ), [ 1, 2 ] ];
        Tjs = Tj[ ( end - 4 ) : ( end ), 4 ];
        ArrVel = [ ];
        for Fr in 2:5
            ( V, ΔT ) = Velocities( Tjs, Locs, Fr );
            push!( ArrVel, V );
        end
        Tjs = Tjs[ 2:end ];
        X2 = zeros( 4, 2 );
        X2[ :, 1 ] = transpose( Tjs );
        X2[ :, 2 ] .= 1.0;
        ArrVel = Float64.( ArrVel );
        coeff_pred = X2 \ ArrVel;
        tol_dist = ( coeff_pred[ 1 ] * ( Tjs[ end ] + 1 ) + coeff_pred[ 2 ] ) * 1.5;
        if tol_dist < ( tol_dist / 5 );
            tol_dist = ( tol_dist / 5 );
        end
    end
    Last = Tj[ end, : ];
    dist = DistanceVectors( Last, cms[ tiempo + 1 ] );
    Closer = argmin( dist );
    DistCloser = minimum( dist );
    if DistCloser < tol_dist
        Temporal = [ transpose( cms[ tiempo + 1 ][ Closer, : ] ) tiempo + 1 ];
        Chain = vcat( Tj, Temporal );
        cms[ tiempo + 1 ] = cms[ tiempo + 1 ][ 1:end .!= Closer, : ];
        t_aux = 0;
    else
        Chain = Tj;
        if t_aux == tol_time
            More = false;
            if size( Tj )[ 1 ] > Tmin
                allTjs[ nTs ] = Tj;
                nTs += 1;
            end
        end
        t_aux += 1;
    end
    return Chain, cms, More, t_aux, allTjs, nTs
end

"""
    Velocities( Tjs, Locs, Fr ) -> V, ΔT
        using DistanceCoords
"""
function Velocities( Tjs, Locs, Fr )
    x = Locs[ Fr - 1, : ];
    y = Locs[ Fr, : ];
    ΔLoc = DistanceCoords( x, y );
    ΔT = Tjs[ Fr ] - Tjs[ Fr  - 1 ];
    V = ΔLoc / ΔT;
    return V, ΔT
end

"""
    FixingGaps( Tjs, Starts, Stops ) -> Tjs
"""
function FixingGaps( Tjs, Starts, Stops )
    for i in 1:length( Tjs )
        if ( Stops[ i ] - Starts[ i ] ) + 1 != size( Tjs[ i ] )[ 1 ]
            Frs = Tjs[ i ][ :, end ];
            for j in 1:( size( Tjs[ i ] )[ 1 ] - 1 )
                J = j + 1;
                dif = Frs[ J ] - Frs[ j ];
                if dif > 1
                    Head = Tjs[ i ][ 1 : j, : ];
                    Body = Tjs[ i ][ J : end, : ];
                    for k in 1: ( dif - 1 )
                        K = ( k / dif );
                        h1 = ( Tjs[ i ][ j, 1 ] + ( Tjs[ i ][ J, 1 ] - Tjs[ i ][ j, 1 ] ) * K );
                        h2 = ( Tjs[ i ][ j, 2 ] + ( Tjs[ i ][ J, 2 ] - Tjs[ i ][ j, 2 ] ) * K );
                        h3 = ( Tjs[ i ][ j, 3 ] + ( Tjs[ i ][ J, 3 ] - Tjs[ i ][ j, 3 ] ) * K );
                        h4 = ( Starts[ i ] + j + k - 1 );
                        Extra = [ h1 h2 h3 h4 ];
                        Head = vcat( Head, Extra );
                    end
                    Tjs[ i ] = vcat( Head, Body );
                end
            end
        end
    end
    return Tjs
end

"""
    StartStop( Tjs ) -> t0, tN
"""
function StartStop( Tjs )
    t0 = Dict{ Int, Float64 }( );
    for i in 1:length( Tjs )
        t0[ i ] = Tjs[ i ][ 1, end ];
    end
    tN = Dict{ Int, Float64 }( );
    for i in 1:length( Tjs )
        tN[ i ] = Tjs[ i ][ end, end ];
    end
    return t0, tN
end
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
end # module STEP02_v1
# ----------------------------------------------------------------------------------------- #
