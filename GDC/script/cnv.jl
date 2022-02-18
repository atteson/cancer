using GDC
using Dates

filters = (Field("files.access") == "open") &
    (Field("files.data_category") == "copy number variation") &
    (Field("files.experimental_strategy") == "WGS")

result = find_files( filters, size=10_000, format="tsv" );

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
            if !isfile( joinpath( cancerdir, "cnv_raw", result[i,:file_id], filename ) )
                push!( files, result[i,:file_id] )
            end
        end
        i += 1
    end
    get_tgz_files( files )
end

@assert( length(setdiff( result[!,:file_id], readdir(joinpath(cancerdir,"cnv_raw")) )) == 0 )

