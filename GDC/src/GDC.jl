module GDC

export Field, GDCString

using HTTP
using CSV
using DataFrames

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

function GDCString( e::Expression, fields::Vector{String}; format::String="csv", size=Inf )
    lines = GDCLines( e )
    string = "{\n\t\"filters\":" * join( lines, "\n\t" ) * ",\n"
    string *= "\t\"format\":\"$format\",\n"
    string *= "\t\"fields\":\"$(join(fields,","))\""
    if size < Inf
        string *= ",\n\t\"size\":\"$size\""
    end
    string *= "\n}"
    return string
end

end


                           
