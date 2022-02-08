using GDC
using Commas
using Statistics

@time snv = read( joinpath( cancerdir, "snv" ), Comma );

casetype = eltype(snv[:case_id])
casecounts = Dict{casetype, Int}()
[casecounts[case] = get( casecounts, case, 0 ) + 1 for case in snv[:case_id]]

findmax(collect(values(casecounts)))

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
