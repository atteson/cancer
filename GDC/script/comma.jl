using GDC
using CSV
using DataFrames
using Commas
using Dates
using StringViews

const newline = UInt8('\n')
const tab = UInt8('\t')
const pound = UInt8('#')

reset( ::Type{Nothing} ) = nothing
reset( ::Type{Int} ) = 0
reset( ::Type{String} ) = UInt8[]

function newlines( buffer, start )
    n = 0
    for i = start:length(buffer)
        n += buffer[i] == newline
    end
    n += buffer[end] != newline
    return n
end

function readmaf( filename )
    numstates = 6

    transitions = zeros( UInt8, numstates, 256 )
    transitions[1,:] = fill( 1, 256 )
    transitions[1,pound] = 2
    transitions[1,tab] = 3
    transitions[1,newline] = 6

    transitions[2,:] = fill( 2, 256 )
    transitions[2,newline] = 1

    transitions[3,:] = fill( 1, 256 )

    transitions[4,:] = fill( 4, 256 )
    transitions[4,tab] = 5
    transitions[4,newline] = 6

    transitions[5,:] = fill( 4, 256 )
    transitions[5,tab] = 5
    transitions[5,newline] = 6
    
    transitions[6,:] = fill( 4, 256 )
    transitions[6,tab] = 5
    transitions[6,newline] = 6

    row = 0
    col = 1
    headers = String[]
    header = UInt8[]
    
    state = 1

    locations = Matrix{Int}( undef, (0, 0) )

    buffer = read( filename )
    for i = 1:length(buffer)
        state = transitions[state,buffer[i]]
        if state == 1
            push!( header, buffer[i] )
        elseif state == 3
            push!( headers, String(header) )
            header = UInt8[]
        elseif state == 5
            col += 1
            locations[row,col] = i
        elseif state == 6
            if row == 0
                push!( headers, String(header) )
                locations = zeros( Int, newlines(buffer, i+1), length(headers)+1 )
            else
                locations[row,col+1] = i
            end
            row += 1
            col = 1
            if row <= size(locations,1)
                locations[row,col] = i
            end
        end
    end
    if row == size(locations,1)
        locations[row,col+1] = i
    else
        @assert( row == size(locations,1)+1 )
    end
    return (buffer, headers, locations)
end

coltypes = Dict(
    "Start_Position" => Int,
    "End_Position" => Int,
    "Start_Position" => Int,
    "End_Position" => Int,
    "t_depth" => Int,
    "t_ref_count" => Int,
    "t_alt_count" => Int,
)

function parseall( buffer, starts, stops, ::Type{String} )
    ranges = UnitRange.( starts .+ 1, stops .- 1 );
    views = view.( [buffer], ranges );
    results = StringView.( views );
    return results
end

function parseall( buffer, starts, stops, ::Type{T} ) where {T}
    results = Vector{T}( undef, length(starts) )
    for i = 1:length(starts)
        range = starts[i] + 1:stops[i] - 1;
        results[i] = parse( T, StringView( view( buffer, starts[i]+1:stops[i]-1 ) ) )
    end
    return results
end

function DataFrames.DataFrame( coltypes, buffer, headers, locations )
    df = DataFrame()
    for i = 1:length(headers)
        header = headers[i]
        coltype = get( coltypes, header, String )
        starts = view( locations, :, i )
        stops = view( locations, :, i+1 )
        df[!,header] = parseall( buffer, starts, stops, coltype );
    end
    return df
end

function get_all( ; printevery = Second(1), dfs = Dict{String,DataFrame}(), maxdir = Inf )
    t0 = now()
    println( "Starting at $t0" )

    rawdir = joinpath( cancerdir, "raw" )
    dirs = readdir( rawdir );

    maxdir = Int(min( length(dirs), maxdir ))
    for i = 1:maxdir
        fileid = dirs[i]

        if haskey( dfs, fileid )
            continue
        end
        
        dir = joinpath( rawdir, fileid )
        filenames = readdir(dir)
        @assert( length(filenames) == 1 )
        filename = joinpath( dir, filenames[1] )

        (buffer, headers, locations) = readmaf( filename )
        dfs[fileid] = df = DataFrame( coltypes, buffer, headers, locations );
        df[!,:file_id] = fill( fileid, size(dfs[fileid],1) )

        m = match( r"^TCGA\.[^\.]*\.([^\.]*)\..", filenames[1] )
        if m != nothing
            caller = m.captures[1]
            if "caller" in names(df)
                @assert( df[!,:callers] .== caller )
            else
                df[!,:callers] = fill( caller, size(df,1) )
            end
        end

        t1 = now()
        if t1 - t0 >= printevery
            println( "Finished $i at $t1" )
            t0 = t1
        end
    end
    return dfs
end

function annotate!( dfs; printevery = Second(1) )
    t0 = now()
    println( "Starting at $t0" )
    
    caseids = union(getindex.( values(dfs), !, :case_id )...);

    increment = 100
    start = 0

    hitdicts = Dict{String,Dict{String,Any}}[]

    while start < length(caseids)
        range = start+1:min(start+increment,length(caseids))
        newcases = endpoint( "cases", Field("case_id") in caseids[range], size = increment+1 )
        hits = newcases["data"]["hits"];
    
        push!( hitdicts, Dict( [hit["case_id"] => hit for hit in hits] ) )
        
        start += increment
        
        t1 = now()
        if t1 - t0 >= printevery
            println( "Finished $start at $t1" )
            t0 = t1
        end
    end
    hitdict = merge(hitdicts...);
    [df[!,:primary_site] = Vector{String}( getindex.( getindex.( [hitdict], df[!,:case_id] ), "primary_site" ) ) for df in values(dfs)];
    [df[!,:disease_type] = Vector{String}( getindex.( getindex.( [hitdict], df[!,:case_id] ), "disease_type" ) ) for df in values(dfs)];
    return hitdict
end

fillers = Dict(
    String1 => String1(" "),
    String31 => String31( "" ),
    String => "",
    Float64 => NaN,
)

unmissing( a::Vector{Union{Missing,T}}, miss::Vector{Bool} ) where {T} =
    Vector{T}( ifelse.( miss, fillers[T], a ) )

function combine( dfs::Dict{String,DataFrame} )
    ns = intersect(names.(values(dfs))...);
    nsdfs = getindex.( values(dfs), !, [ns] );

    df = DataFrame()
    for i = 1:length(ns)
        println( "Processing $(ns[i]) at $(now())" )
        as = getindex.( nsdfs, !, ns[i] );
        a = vcat( as... );
        miss = ismissing.(a);
        if sum(miss) > 0
            @time a = unmissing(a, miss);
        end
        df[!,ns[i]] = a
    end
    return df
end

function writecomma( df )
    newdf = DataFrame()
    N = size(df,1)
    for n in names(df)
        println( "Processing $n at $(now())" )
        if eltype(df[!,n]) <: AbstractString
            l = mapreduce( length, max, df[!,n] )
            if l > 0
                v = fill( UInt8(' '), N*l )

                j = 1
                for i = 1:N
                    copyto!( v, j, df[i,n] )
                    j += l
                end
                newdf[!,n] = reinterpret( CharN{l}, v )
            end
        else
            newdf[!,n] = df[!,n]
        end
    end
    return newdf
end

@time dfs = get_all();

@time hitdict = annotate!( dfs );
                            
@time df = combine( dfs );

@time comma = Comma( joinpath( cancerdir, "snv_new" ), df );



@time df1 = tocharn( df );

@time c = Comma( df1 );

@time write( joinpath( cancerdir, "snv" ), c )

@time c2 = read( joinpath( cancerdir, "snv" ), Comma )

c = read( joinpath( cancerdir, "snv_new" ), Comma );
