using Distributions

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
    return probs[n % 1 + 1,:]
end

@time p10 = calculate( 0.9, 10 );
@time p12 = calculate( 0.9, 12 );
@time p15 = calculate( 0.9, 15 );
