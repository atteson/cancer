using HTTP
using CSV
using DataFrames

struct Field
    name::Symbol
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

e = ( Field(:x) == 2 ) & ( Field(:y) == 3 ) & ( Field(:z) in [:a,:b,:c] )

attributes = Dict( [
    :(==) => ["=","{","}"],
    :(&) => ["and","[","]"],
    :in => ["in","[","]"],
    :(|) => ["or","[","]"],
] )

GDCString( x ) = join( GDCLines( x ), "\n" )

function GDCLines( e::Expression{op} ) where {op}
    (opstring, start, stop) = attributes[op]
    prefix = [ "{",
               "\t\"op\":\"$(transform[op])\"",
               "\t\"content\":$start" ]
                
    body = GDCLines.( e.args )
    body = [vcat(map( v -> [v[1:end-1]; v[end] * ","], body[1:end-1] )...); body[end]]
               
    suffix = [ "\t$stop", "}" ]
    return [prefix; "\t\t" .* body; suffix]
end
    
GDCLines( field::Field ) = ["\"field\":\"$(field.name)\""]

GDCLines( v::Vector ) = [ "\"value\":["; "\t\"" .* [v[1:end-1] .* "\","; v[end] * "\""]; "]" ]

GDCLines( x ) = ["\"value\":\"$x\""]

println(GDCString( (Field(:x) == 2) & (Field(:y) == 3) & (Field(:z) in ["a","b","c"]) ))

                           
