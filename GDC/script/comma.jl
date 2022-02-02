using GDC
using CSV
using DataFrames
using Commas
using Dates

fillers = Dict(
    String1 => String1(" "),
    String31 => String31( "" ),
    String => "",
    Float64 => NaN,
)

function get_all( ; printevery = Second(1), start=1, stop=Inf, df = DataFrame() )
    cancerdir = joinpath( homedir(), "data", "cancer" )

    dirs = readdir( cancerdir );

    stop = Int(min(length(dirs),stop))
    
    t0 = now()
    println( "Starting at $t0" )
    for i = start:stop-1

        fileid = dirs[i]

        if !isempty(df) && fileid in df[!,:file_id]
            continue
        end
        
        dir = joinpath( cancerdir, fileid )
        filenames = readdir(dir)
        @assert( length(filenames) == 1 )
        filename = joinpath( dir, filenames[1] )
        data = CSV.read( filename, DataFrame, comment="#", delim="\t" );
        
        data[!,:file_id] = fill( fileid, size(data,1) )

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

        t1 = now()
        if t1 - t0 >= printevery
            println( "Finished $i at $t1" )
            t0 = t1
        end
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

        hits = cases["data"]["hits"]
        @assert( sort(getindex.( hits, "case_id" )) == sort(caseids) )

        hitdict = Dict( [hit["case_id"] => hit for hit in hits] )
        casehits = getindex.( [hitdict], data[!,:case_id] );
                   
        data[!,:disease_type] = getindex.( casehits, "disease_type" )
        data[!,:primary_site] = getindex.( casehits, "primary_site" )
end
