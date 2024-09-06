module GDC

export Field, find_files, get_files, endpoint, cancerdir, get_tgz_files

using HTTP
using CSV
using DataFrames
using Tar
using CodecZlib
using JSON

const cancerdir = joinpath( homedir(), "data", "cancer" )

struct Field
    name::String
end

struct Expression{T}
    args::Vector
end

Base.:(==)( f::Field, x ) = Expression{:(==)}( [f, x] )
Base.:(==)( x, f::Field ) = Expression{:(==)}( [f, x] )

# && and || are syntatic since they are short circuit
Base.:(&)( e1::Expression, e2::Expression ) = Expression{:(&)}( [e1, e2] )
Base.:(|)( e1::Expression, e2::Expression ) = Expression{:(|)}( [e1, e2] )

Base.in( f::Field, x ) = Expression{:in}( [f, x] )

attributes = Dict( [
    :(==) => ["=","{","}"],
    :(&) => ["and","[","]"],
    :in => ["in","{","}"],
    :(|) => ["or","[","]"],
] )

function GDCLines( e::Expression{op} ) where {op}
    (opstring, start, stop) = attributes[op]
    prefix = [ "{",
               "\t\"op\":\"$opstring\",",
               "\t\"content\":$start" ]
                
    body = GDCLines.( e.args )
    body = [vcat(map( v -> [v[1:end-1]; v[end] * ","], body[1:end-1] )...); body[end]]
               
    suffix = [ "\t$stop", "}" ]
    return [prefix; "\t\t" .* body; suffix]
end
    
GDCLines( field::Field ) = ["\"field\":\"$(field.name)\""]

# quotes in quotes seems to screw up emacs syntax highlighting
q = '"'
GDCLines( v::Vector ) = [ "$(q)value$q:["; "\t$q" .* [v[1:end-1] .* "$q,"; v[end] * "$q"]; "]" ]

GDCLines( x ) = ["\"value\":\"$x\""]

function GDCString( e::Expression, fields::Vector{String}; kwargs... )
    lines = GDCLines( e )

    # the value corresponding to filters is not a string so we do this one separately
    s = "{\n\t\"filters\":" * join( lines, "\n\t" )*","
    
    kvs = (fields = join( fields, "," ), kwargs... )
    s *= join(["\n\t\"$k\":\"$v\"" for (k,v) in pairs(kvs)], ",")
    s *= "\n}"
    return s
end

function endpoint( endpoint, filter; fields = String[], kwargs... )
    body = GDCString( filter, fields; kwargs... )
    response = HTTP.request( "POST", "https://api.gdc.cancer.gov/$endpoint", ["Content-Type" => "application/json"], body );
    return JSON.parse( String(response.body) )
end

file_fields = [
    "file_id",
    "file_name",
    "data_type",
    "data_format",
]

function find_files( filter::Expression; fields = file_fields, kwargs... )
    body = GDCString( filter, fields; kwargs... )
    response = HTTP.request( "POST", "https://api.gdc.cancer.gov/files", ["Content-Type" => "application/json"], body );
    return CSV.read( response.body, DataFrame );
end

function request_files( files::Vector{String}; debug=false )
    url = "https://api.gdc.cancer.gov/data"
    body = "{\n\t\"ids\":[\n\t\t\"" * join( files, "\",\n\t\t\"" ) * "\"\n\t]\n}"
    debug && print( "Posting to $url with body:")
    debug && print( body )
    response = HTTP.request( "POST", url, ["Content-Type" => "application/json"], body );

    filename = match( r"attachment; filename=(.*)$", Dict(response.headers)["Content-Disposition"] ).captures[1]

    return (response, filename)
end
   
function get_tgz_files( files::Vector{String}; debug=false, dir = joinpath( cancerdir, "raw" ) )
    debug && println( "Requesting files" )
    (response, filename) = request_files( files, debug=debug )
    debug && println( "Done requesting files" )

    debug && println( "Extracting files" )
    tmpdir = tempname()
    stream = GzipDecompressorStream( IOBuffer(response.body) )
    if occursin( r"\.tar\.gz$", filename )
        Tar.extract( stream, tmpdir )
    else
        path = joinpath( tmpdir, split( filename, "." )[1] )
        mkpath( path )
        fileroot = match( r"^(.*)\.gz$", filename ).captures[1]
        write( joinpath( path, fileroot ), stream )
    end
    debug && println( "Done extracting files" )
    debug && println( "Copying files from temp" )
    files = setdiff( readdir(tmpdir), ["MANIFEST.txt"] )
    cp.( joinpath.( tmpdir, files ), joinpath.( dir, files ), force=true )
    rm( tmpdir, recursive=true )
    debug && println( "Done copying files from temp" )
end

function uncompress( filename::AbstractString )
end

end


                           
