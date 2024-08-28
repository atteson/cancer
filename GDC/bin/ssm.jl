using GDC
using Commas
using StatsBase
using Statistics
using DataFrames
using CSV
using Plots
using DataStructures

@time snv = read( joinpath( cancerdir, "snv" ), Comma );
@time cnv = read( joinpath( cancerdir, "cnv" ), Comma );

snv_cases = unique(snv[:case_id]);
cnv_cases = unique(cnv[:case_id]);

cases = snv_cases[in.(snv_cases,[cnv_cases])];
length(cases)

@time good = snv[in.(snv[:case_id],[Set(cases)])];

scnv = sort( cnv, :case_id, :Chromosome, :Start );

mean(cnv[:Copy_Number] .== cnv[:Major_Copy_Number])

sgood = sort( good, :case_id, :Chromosome, :Start_Position );

function find_indices( sgood, scnv )
    i = 1
    j = 1
    indices = Int[]
    sgoodstart = getindex.( [sgood], i, (:case_id,:Chromosome,:Start_Position) )
    scnvstart = getindex.( [scnv], j, (:case_id,:Chromosome,:Start) )
    scnvend = getindex.( [scnv], j, (:case_id,:Chromosome,:End) )
    while i <= size(sgood,1) && j <= size(scnv,1)
        if sgoodstart >= scnvstart
            if sgoodstart <= scnvend
                push!( indices, j )
                i += 1
                if i <= size(sgood,1)
                    sgoodstart = getindex.( [sgood], i, (:case_id,:Chromosome,:Start_Position) )
                end
            else
                j += 1
                if j <= size(scnv,1)
                    scnvstart = getindex.( [scnv], j, (:case_id,:Chromosome,:Start) )
                    scnvend = getindex.( [scnv], j, (:case_id,:Chromosome,:End) )
                end
            end
        else
            push!( indices, 0 )
            i += 1
            if i <= size(sgood,1)
                sgoodstart = getindex.( [sgood], i, (:case_id,:Chromosome,:Start_Position) )
            end
        end
    end
    while i <= size(sgood,1)
        push!( indices, 0 )
        i += 1
    end
    return indices
end

@time indices = find_indices( sgood, scnv );

@assert( length(indices) == size(sgood,1) )

sum(indices.==0)

size(sgood)

@time cases = [sgood; :index => indices];
@time cases = cases[cases[:index] .!= 0];
@time cases = [cases;
               :Minor_Copy_Number => scnv[cases[:index],:Minor_Copy_Number];
               :Major_Copy_Number => scnv[cases[:index],:Major_Copy_Number]
               ];

loss = (cases[:Minor_Copy_Number].==0) .& (cases[:Major_Copy_Number].==1);
chrX = convert( CharN{5}, "chrX" )
x = filter( x->x[1] == chrX,
            sort(collect(zip(cases[loss,:Chromosome], cases[loss,:t_alt_count]./cases[loss,:t_depth])))
            )
histogram( getindex.( x, 2 ) )

case = collect( groupby( cases, :case_id ) )[1];

@assert( length(unique(zip(case[:Chromosome], case[:Start_Position]))) == size(case,1) )
@assert( all(case[:index] .!= 0) )

rundir = joinpath( homedir(), "cancer", "pairtree", "runs" )
mkpath( rundir )
dir = tempname( rundir )
mkpath( dir )
file = open( joinpath( dir, "case.ssm" ), "w" )
write( file, join( [:id,:name,:var_reads,:total_reads,:var_read_prob], "\t" ) )

ids = "s" .* string.(0:size(case,1)-1)
names = convert.( String, case[:Chromosome] ) .* "_" .* string.(case[:Start_Position])
var_reads = case[:t_alt_count]
total_reads = case[:t_alt_count] .+ case[:t_ref_count]
var_read_probs
collect(zip( scnv[case[:index],:Minor_Copy_Number], scnv[case[:index],:Major_Copy_Number],
             case[:t_alt_count]./case[:t_depth], case[:Chromosome]))

case[:t_alt_count] .+ case[:t_ref_count] .- case[:t_depth]
    
close( file )


comma = case;
@time c = Comma(NamedTuple{keys(comma)}(collect.( getindex.( [comma], keys(comma) ) )));


