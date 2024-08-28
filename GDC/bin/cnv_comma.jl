using GDC
using CSV
using DataFrames
using Commas
using Dates
using StringViews

function get_all( ; printevery = Second(1), dfs = Dict{String,DataFrame}(), maxdir = Inf )
    t0 = now()
    println( "Starting at $t0" )

    rawdir = joinpath( cancerdir, "cnv_raw" )
    mkpath( rawdir )
    dirs = readdir( rawdir );

    maxdir = Int(min( length(dirs), maxdir ))
    for i = 1:maxdir
        fileid = dirs[i]

        dir = joinpath( rawdir, fileid )
        if haskey( dfs, fileid ) || !isdir( dir )
            continue
        end
        
        filenames = readdir(dir)
        @assert( length(filenames) == 1 )
        filename = joinpath( dir, filenames[1] )
        if !occursin( r"\.wgs\.ASCAT\.copy_number_variation\.seg\.txt$", filename )
            continue
        end

        df = dfs[fileid] = CSV.read( filename, DataFrame, delim='\t', comment="#", pool=false );
        df[!,:file_id] = fill( fileid, size(dfs[fileid],1) )

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
    
    fileids = union(getindex.( values(dfs), !, :file_id )...);
    aliquotids = union(getindex.( values(dfs), !, :GDC_Aliquot )...);

    increment = 100
    start = 0

    aliquotdicts = Dict{String,String}[]

    while start < length(fileids)
        range = start+1:min(start+increment,length(fileids))

        newcases = endpoint( "cases", Field("files.file_id") in fileids[range], size=increment+1 )
        
        hits = newcases["data"]["hits"];

        newaliquotids = getindex.( hits, "aliquot_ids" )
        newcaseids = fill.( getindex.( hits, "case_id" ), length.(newaliquotids) )
        push!( aliquotdicts, Dict( union( [newaliquotids[i] .=> newcaseids[i] for i in 1:length(newcaseids)]... )... ) )

        @assert( isempty(setdiff( aliquotids[range], keys(aliquotdicts[end]) )) )
        
        start += increment
        
        t1 = now()
        if t1 - t0 >= printevery
            println( "Finished $start at $t1" )
            t0 = t1
        end
    end

    aliquotdict = merge(aliquotdicts...);
    return aliquotdict
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

@time dfs = get_all();

@time aliquotdict = annotate!( dfs );
                            
@time df = combine( dfs );

df[!,:case_id] = getindex.( [aliquotdict], df[!,:GDC_Aliquot] )

@time comma = Comma( joinpath( cancerdir, "cnv_new" ), df );

c = read( joinpath( cancerdir, "cnv_new" ), Comma );
