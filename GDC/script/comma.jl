using GDC
using CSV
using DataFrames
using Commas
using Dates

function get_all( ; printevery = Second(1), start=1, stop=Inf, dfs = Dict{String,DataFrame}() )
    cancerdir = joinpath( homedir(), "data", "cancer" )
    
    dirs = readdir( cancerdir );

    stop = Int(min(length(dirs),stop))
    
    t0 = now()
    println( "Starting at $t0" )
    for i = start:stop
        fileid = dirs[i]

        if haskey( dfs, fileid )
            continue
        end
        
        dir = joinpath( cancerdir, fileid )
        filenames = readdir(dir)
        @assert( length(filenames) == 1 )
        filename = joinpath( dir, filenames[1] )
        dfs[fileid] = CSV.read( filename, DataFrame, comment="#", delim="\t" );
        dfs[fileid][!,:file_id] = fill( fileid, size(dfs[fileid],1) )

        t1 = now()
        if t1 - t0 >= printevery
            println( "Finished $i at $t1" )
            t0 = t1
        end
    end
    return dfs
end

function annotate!( dfs; printevery = Second(1) )
    caseids = union(getindex.( values(dfs), !, :case_id )...);

    increment = 100
    start = 0

    hitdicts = Dict{String,Dict{String,Any}}[]

    t0 = now()
    println( "Starting at $t0" )
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
    [df[!,:primary_site] = getindex.( getindex.( [hitdict], df[!,:case_id] ), "primary_site" ) for df in values(dfs)];
    [df[!,:disease_type] = getindex.( getindex.( [hitdict], df[!,:case_id] ), "disease_type" ) for df in values(dfs)];
    return hitdict
end

dfs = get_all();

hitdict = annotate!( dfs );

fillers = Dict(
    String1 => String1(" "),
    String31 => String31( "" ),
    String => "",
    Float64 => NaN,
)

unmissing( a::Vector{Union{Missing,T}}, miss::Vector{Bool} ) where {T} =
    Vector{T}( ifelse.( miss, fillers[T], a ) )

function combine( dfs )
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
end

function f()
        df2 = data
        if !isempty(df)
            fields = names(df)
            ns = names(data)
            extra = length(setdiff(ns, fields))
            if extra > 0
                println( "Data set $i has $extra extra columns" )
            end

            miss = setdiff(fields, ns)
            if !isempty(miss)
                println( "Data set $i is missing $(length(miss)) columns" )

                types = [unique(typeof.(df[!,n])) for n in miss]
                types = [t[t .!= Missing] for t in types]
                @assert( unique(length.(types)) == [1] )

                types = getindex.( types, 1 )

                setindex!.( [data], fill.( getindex.( [fillers], types ), size(data,1) ), !, miss )
            end
            df2 = data[!,fields]
        end

        df = vcat( df, df2 )
    end
    return df
end

df = DataFrame()
df = get_all( df=df );

function get_cases( df )
        caseids = unique(data[!,:case_id])

        cases = try
            endpoint( "cases", Field("case_id") in caseids, size = 2*length(caseids) )
        catch e
            println( e )
            return df
        end


        @assert( sort(getindex.( hits, "case_id" )) == sort(caseids) )

        hitdict = Dict( [hit["case_id"] => hit for hit in hits] )
        casehits = getindex.( [hitdict], data[!,:case_id] );
                   
        data[!,:disease_type] = getindex.( casehits, "disease_type" )
        data[!,:primary_site] = getindex.( casehits, "primary_site" )
end
