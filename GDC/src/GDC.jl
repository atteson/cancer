module GDC

export Field, find_files, get_files

using HTTP
using CSV
using DataFrames
using Tar
using CodecZlib

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

GDCLines( v::Vector ) = [ "\"value\":["; "\t\"" .* [v[1:end-1] .* "\","; v[end] * "\""]; "]" ]

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

file_fields = [
    "file_id",
    "file_name",
    "data_type",
    "data_format",
]

function find_files( filter::Expression; kwargs... )
    body = GDCString( filter, file_fields; kwargs... )
    response = HTTP.request( "POST", "https://api.gdc.cancer.gov/files", ["Content-Type" => "application/json"], body );
    return CSV.read( response.body, DataFrame );
end

function get_files( files::Vector{String} )
    body = "{\n\t\"ids\":[\n\t\t\"" * join( files, "\",\n\t\t\"" ) * "\"\n\t]\n}"
    response = HTTP.request( "POST", "https://api.gdc.cancer.gov/data", ["Content-Type" => "application/json"], body );

    filename = match( r"attachment; filename=(.*)$", Dict(response.headers)["Content-Disposition"] ).captures[1]
    
    dir = tempname()
    stream = GzipDecompressorStream( IOBuffer(response.body) )
    if occursin( r"\.tar\.gz$", filename )
        Tar.extract( stream, dir )
    else
        path = joinpath( dir, split( filename, "." )[1] )
        mkpath( path )
        fileroot = match( r"^(.*)\.gz$", filename ).captures[1]
        write( joinpath( path, fileroot ), stream )
    end
    files = setdiff( readdir(dir), ["MANIFEST.txt"] )
    cp.( joinpath.( dir, files ), joinpath.( homedir(), "data", "cancer", files ) )
    rm( dir, recursive=true )
end

end


                           
