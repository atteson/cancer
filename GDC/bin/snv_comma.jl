using GDC
using CSV
using DataFrames
using Commas
using Dates
using StringViews
using GZip
using SentinelArrays
using MissingTypes

coltypes = Dict(
    "Start_Position" => Int,
    "End_Position" => Int,
    "Start_Position" => Int,
    "End_Position" => Int,
    "t_depth" => Int,
    "t_ref_count" => Int,
    "t_alt_count" => Int,
    "Tumor_Seq_Allele2" => String,
    "Hugo_Symbol" => String,
)

variable_length_string_cols = ["all_effects"]

function get_all( ; printevery = Second(1), dfs = Dict{String,DataFrame}(), maxdir = Inf )
    t0 = now()
    println( "Starting at $t0" )

    rawdir = joinpath( cancerdir, "snv_raw" )
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

        dfs[fileid] = df = CSV.read( filename, DataFrame, pool=false, delim='\t', comment="#",
                                     types = coltypes );
        
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

function combine( dfs::Dict{String,DataFrame} )
    ns = intersect(names.(values(dfs))...);
    nsdfs = getindex.( values(dfs), !, [ns] );

    df = DataFrame()
    for i = 1:length(ns)
        println( "Processing $(ns[i]) at $(now())" )
        as = getindex.( nsdfs, !, ns[i] );
#        as = as[length.(as).!=0]
#        as = promote(as...)
        a = vcat( as... );
        df[!,ns[i]] = a
    end
    return df
end

@time dfs = get_all( maxdir = 1_000 );

@time hitdict = annotate!( dfs );
                            
@time df = combine( dfs );

dir = joinpath( cancerdir, "snv_new" )
rm( dir, recursive=true )
@time  Comma( dir, df, verbose=true );

read( dir, Comma )

@time ss = [string.(df[!,:all_effects]) for df in values(dfs)];
@time s = vcat(ss...);
@time n = sum(length.(s));
@time v = Vector{UInt8}( undef, n );

function c( vs::Vector{String}, vu::Vector{UInt8} )
    i = 1
    for j = 1:length(vs)
        s = vs[j]
        n = length(s)
        k = 1
        while k <= n
            vu[i] = s[k]
            i += 1
            k += 1
        end
    end
end

@time c( s, v );

@time s = reduce( *, s );
@time s = reduce( *, reduce.( *, ss ) );
