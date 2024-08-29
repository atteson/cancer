using Distributions
using Dates
using SparseArrays

function calculate( p, n )
    N = 2^n

    odds = p/(1-p)
    probs = zeros(2, N)
    probs[1,1] = 1.0
    @assert(sum(probs, dims=2)== [1.0 0.0]')

    bins = Binomial.( 1:2^(n-1), p );

    currsize = 1
    for i = 1:n
        next = i % 2 + 1
        curr = 3 - next

        nextsize = 2 * currsize
        for j = 1:nextsize
            probs[next, j] = 0.0
            for k = Int(ceil( j/2 )):min(j, currsize)
                probs[next, j] += pdf( bins[k], j - k ) * probs[curr, k]
            end
        end
        currsize = nextsize
    end
    return probs[n % 2 + 1,:]
end

function calculate2( p, n )
    No2 = 2^(n-1)
    M = spzeros(No2 << 1, No2)

    bins = Binomial.( 1:No2, p )

    for from = 1:No2
        for to = from:from << 1
            M[to, from] = pdf( bins[from], to - from )
        end
    end

    return M
    
    ps = Vector{Float64}[]
    p = ones(1)
    currsize = 1
    for i = 1:n
        nextsize = currsize << 1
        p = view( M, 1:nextsize, 1:currsize ) * p
        push!( ps, p )
    end
    
    return ps
end

ps = Vector{Float64}[]
for i = 1:16
    println( "Running $i at $(now())" )
    push!( ps, calculate( 0.9, i ) )
end

using Plots

function plot_densities( ps, r; factor=1.0, p = plot( size=[1000,800] ), xscale=1.0, yscale=1.0 )
    for i in r
        pn = ps[i]
        n = length(pn)
        dx = 1/(factor ^ i  * n)
        plot!(p, xscale.*(1:n).*dx, yscale.*pn./dx, label=string(i))
    end
    display(p)
    return p
end

plot_densities( ps, 11:16 )
plot_densities( ps, 11:16, factor=0.95 )
plot_densities( ps, 11:16, factor=sqrt(0.9) )

p = plot_densities( ps, 16:16 )
p = plot_densities( ps, 16:16, xscale=0.5, p=p, yscale=0.25 )

M = calculate2( 0.9, 2 )







