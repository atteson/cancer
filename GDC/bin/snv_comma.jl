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
    "Allele" => String,
    "End_Position" => Int,
    "HGVSp" => String,
    "Hugo_Symbol" => String,
    "PUBMED" => String,
    "PHENO" => String,
    "Reference_Allele" => String,
    "SOMATIC" => String,
    "Start_Position" => Int,
    "Tumor_Seq_Allele1" => String,
    "Tumor_Seq_Allele2" => String,
    "t_alt_count" => Int,
    "t_depth" => Int,
    "t_ref_count" => Int,
)

function get_all( fileids::Set{String}; printevery = Second(1), dfs = Dict{String,DataFrame}(), maxdir = Inf )
    t0 = now()
    println( "Starting at $t0" )

    rawdir = joinpath( cancerdir, "snv_raw" )
    dirs = readdir( rawdir );

    i = 1
    while length(dfs) < maxdir && i <= length(dirs)
        fileid = dirs[i]

        if in( fileid, fileids )
            i += 1
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
        i += 1
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
    for name in ns
        println( "Processing $name at $(now())" )
        as = getindex.( nsdfs, !, name );
        a = vcat( as... );
        df[!,name] = a
    end
    return df
end

dir = joinpath( cancerdir, "snv_new" )
if isdir(dir)
    rm( dir, recursive=true )
end
fileids = Set{String}()

n = 1000
i = 1
while n == 1000
    println( "Processing dataframe $i at $(now())" )
    @time dfs = get_all( fileids, maxdir = n );
    n = length(dfs)

    @time hitdict = annotate!( dfs );
                            
    @time df = combine( dfs );

    @time comma = Comma( df );
    @time write( dir, comma, append=true, verbose=true )

    fileids = union( fileids, string.(unique(comma[:file_id])));
    i += 1
end
