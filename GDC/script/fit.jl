using GDC
using Commas
using StatsBase
using Statistics
using DataFrames
using CSV
using Plots
using DataStructures

@time snv = read( joinpath( cancerdir, "snv" ), Comma );

casetype = eltype(snv[:case_id])
casecounts = Dict{casetype, Int}()
[casecounts[case] = get( casecounts, case, 0 ) + 1 for case in snv[:case_id]];

(m,i) = findmax(collect(values(casecounts)))
badcase = collect(keys(casecounts))[i]

@time df = snv[snv[:case_id].==[badcase]];
size(df)
callertype = eltype(snv[:callers])
callercounts = Dict{callertype, Int}()
[callercounts[caller] = get( callercounts, caller, 0 ) + 1 for caller in df[:callers]];
callercounts

unique(df[:Tumor_Sample_Barcode])

unique(df[:Tumor_Sample_UUID])
startperm = sortperm(df[:Start_Position]);
badstarts = df[:Start_Position][startperm];
badchromes = df[:Chromosome][startperm];
sum((badstarts[1:end-1] .== badstarts[2:end]) .& (badchromes[1:end-1].==badchromes[2:end]))


function multicall( df )
    perm = sortperm( df[:Start_Position] );
    counts = zeros( Int, 2, 2, 2, 2 )

    caller = Dict( zip( unique( df[:callers] ), 1:4 ) )
    index = ones( Int, 4 )
    
    start = 1
    startpos = df[:Start_Position]
    chrome = df[:Chromosome]
    callers = df[:callers]
    index[caller[callers[perm[1]]]] = 2
    for i = 2:length(perm)
        if startpos[perm[i]] == startpos[perm[start]] && chrome[perm[i]] == chrome[perm[start]]
            distinctcallers = true
            for j = start:i-1
                distinctcallers = distinctcallers && callers[perm[j]] != callers[perm[i]]
            end
            if !distinctcallers
                error( "Indistinct callers at $i" )
            end
        else
            counts[index...] = counts[index...] + 1
            start = i
            index = ones( Int, 4 )
        end
        index[caller[callers[perm[i]]]] = 2
    end
    counts[index...] = counts[index...] + 1
    return (caller, counts)
end

(caller, counts) = multicall( df )

cross( x, y ) = vec([[x[i]; y[j]] for i in 1:length(x), j in 1:length(y)])
S = [[1], [2]]
indices = reduce( cross, fill( S, 4 ) )
@assert( sum([counts[index...]*sum(index .- 1) for index in indices]) == size(df,1) )

rellac = strip.(string.(getindex.( [Dict( zip( values(caller), keys(caller) ) )], 1:4 )))
result = DataFrame()
[result[!,caller] = Bool[] for caller in rellac]
result[!,"number of callers"] = Int[]
result[!,:count] = Int[]

[push!( result, [index .== 2; sum(index.==2); counts[index...]] ) for index in indices]
CSV.write( joinpath( homedir(), "table.csv" ), result )

T = eltype(df[:dbSNP_RS])
df2 = df[df[:dbSNP_RS] .== [convert( T, "novel" )],:]

(caller, counts) = multicall( df2 )


s = sort(df[:t_alt_count]./(df[:t_ref_count] .+ df[:t_alt_count]))
quantile( s, 0.25 )

histogram( s )

s = sort(df[:t_ref_count] .+ df[:t_alt_count])
quantile( s, 0.5 )
quantile( df[:t_depth], 0.5 )

sdf = sort( df, :Start_Position, :Chromosome );

vafs = [group[:t_alt_count]./(group[:t_alt_count] .+ group[:t_ref_count]) for group in Commas.groupby(sdf) if size(group,1)>1]
sort(std.(vafs))
p = histogram( mean.(vafs), label="" )
savefig( p, joinpath( homedir(), "vafs.png" ) )

ps = sort( snv, :primary_site )
psd = Dict([group[1,:primary_site] => length(unique(group[:case_id])) for group in Commas.groupby(ps)])
perm = sortperm(collect(values(psd)));
collect(zip(collect(keys(psd))[perm], collect(values(psd))[perm]))


cc = sort(collect(values(casecounts)))
quantile( cc, 0.5 )
mean(cc)

unique(snv[:primary_site])
casesites = Dict( [snv[:case_id][i]=>snv[:primary_site][i] for i in 1:size(snv,1)]... );

sitetype = eltype(snv[:primary_site])
sitecounts = Dict{sitetype,Vector{Int}}()
for (case,site) in casesites
    if !haskey( sitecounts, site )
        sitecounts[site] = Int[]
    end
    push!( sitecounts[site], casecounts[case] )
end    

[site => (quantile(counts,0.5),mean(counts),maximum(counts)) for (site, counts) in sitecounts]
