using GDC
using Dates
using Debugger

filters = (Field("files.access") == "open") &
    (Field("files.data_category") == "sequencing reads") &
    (Field("files.experimental_strategy") == "WGS")

result = find_files( filters, size=10_000, format="tsv" );

dir = joinpath( cancerdir, "raw" )

i = 1
increment = 10
n = size(result,1)
while i <= n
    println( "Processing next batch at $(now())" )
    files = String[]
    while length(files) < increment && i <= n
        m = match( r"^(.*\.wgs\.ASCAT\.(gene_level\.copy_number_variation\.tsv|copy_number_variation\.seg\.txt))$",
                   result[i,:file_name] )
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
    sleep(10)
    get_tgz_files( files, debug=true )
end

@assert( length(setdiff( result[!,:file_id], readdir(dir) )) == 0 )
