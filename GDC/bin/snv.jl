using GDC
using Dates

# no WGS for SNV
filters = (Field("files.access") == "open") &
    (Field("files.data_category") == "simple nucleotide variation") &
    (Field("files.experimental_strategy") == "WXS")

result = find_files( filters, size=20_000, format="tsv" );

dir = joinpath( cancerdir, "snv_raw" )
mkpath(dir)

i = 1
increment = 10
n = size(result,1)
while i <= n
    println( "Processing next batch at $(now())" )
    files = String[]
    while length(files) < increment && i <= n
        m = match( r"^(.*\.maf)(\.tar)?(\.gz)?$", result[i,:file_name] )
        if m == nothing
            error( "$(result[i,:file_name]) doesn't match" )
        else
            filename = m.captures[1]
            if !isfile( joinpath( dir, result[i,:file_id], filename ) )
                push!( files, result[i,:file_id] )
            end
        end
        i += 1
    end
    get_tgz_files( files, dir=dir, debug=true )
end

@assert( length(setdiff( result[!,:file_id], readdir(dir) )) == 0 )
