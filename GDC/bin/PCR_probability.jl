using Plots
using Distributions
using Dates
using SparseArrays
using LinearAlgebra

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

function calculate_M( p, n )
    No2 = 2^(n-1)
    M = spzeros(No2 << 1, No2)

    bins = Binomial.( 1:No2, p )

    for from = 1:No2
        for to = from:from << 1
            M[to, from] = pdf( bins[from], to - from )
        end
    end

    return M
end

function calculate_partial_M( p, start, n, maxsize )
    No2 = 1 << (n-1)
    N = No2 << 1
    M = spzeros( N, No2 )
    i = start
    nz = 0
    factor = p/(1-p)

    while nz + i < maxsize && i <= No2
        base = i + Int(round(p*i))
        i2 = i << 1
        M[base, i] = pdf( Binomial( i, p ), base - i )
        
        up = base+1:i2
        down = base-1:-1:i
        upfactor = factor .* (i2 .- up .+ 1)./(up .- i)
        upresult = accumulate( *, upfactor, init=M[base, i] )
        M[up, i] = upresult
        downfactor = (down .+ 1 .- i)./(i2 .- down)./factor
        downresult = accumulate( *, downfactor, init=M[base, i] )
        M[down, i] = downresult

        nz += i2
        i += 1
    end
    return M
end

function calculate_P( p, n; vec=false )
    No2 = 1 << (n-1)
    N = No2 << 1
    i = 1
    factor = p/(1-p)

    P = zeros( n+1, N )
    P[1,1] = 1.0
    
    while i <= No2
        m = zeros( i+1 )
        base = i + Int(round(p*i))
        i2 = i << 1

        baseresult = pdf( Binomial( i, p ), base - i )

        m[base - i + 1] = baseresult
        
        up = base+1:i2
        down = base-1:-1:i
        upfactor = factor .* (i2 .- up .+ 1)./(up .- i)
        upresult = accumulate( *, upfactor, init=baseresult )

        c = length(up)
        if c > 0
            m[up .- i .+ 1] = upresult
        end
        
        downfactor = (down .+ 1 .- i)./(i2 .- down)./factor
        downresult = accumulate( *, downfactor, init=baseresult )

        m[down .- i .+ 1] = downresult

        if !vec
            for k = 1:n
                for l = i:i2
                    P[k+1,l] += P[k,i] * m[l-i+1]
                end
            end
        else
            r = i:i2
            for k = 1:n
                result = P[k+1,r] + P[k,i] * m[r .- i .+ 1]
                P[k+1,r] = result
            end
        end

        i += 1
    end
    return P
end

n = 14
No2 = 1 << (n-1)
N = No2 << 1
P0 = zeros( n+1, N )
P0[1,1] = 1
@time M0 = calculate_M( 0.9, n );
@time for k = 1:n
    P0[k+1,:] = M0 * P0[k,1:No2]
end
@time P1 = calculate_P( 0.9, n );
@time P1 = calculate_P( 0.9, n );
@time P1 = calculate_P( 0.9, n, vec=true );
@time P1 = calculate_P( 0.9, n, vec=true );
maximum(abs.(P1 - P0))

ts = []
for n = 1:14
    t0 = time()
    calculate_P( 0.9, n );
    t1 = time()
    push!( ts, t1 - t0 )
end
plot(ts)
ts[2:end]./ts[1:end-1]
diff(log.(ts))
ts[end]*4^(25-14)

using Profile
using ProfileView
Profile.clear()
@profile P1 = calculate_P( 0.9, n );
ProfileView.view()


@time M1 = calculate_partial_M( 0.9, 1, n );
maximum(abs.(M0 - M1))

using Profile
using ProfileView
Profile.clear()
@profile M1 = calculate_partial_M( 0.9, 1, 11 );
ProfileView.view()

ps = Vector{Float64}[]
for i = 1:16
    println( "Running $i at $(now())" )
    push!( ps, calculate( 0.9, i ) )
end

function plot_densities( ps, r; factor=1.0, p = plot( size=[1000,800] ), xs=1.0, ys=1.0, kwargs... )
    for i in r
        pn = ps[i]
        j = (f -> f == nothing ? length(pn) : f-1 )(findfirst( pn .== 0 ))
        pn = pn[1:j]
        n = length(pn)
        dx = 1/((2 * factor) ^ i )
        plot!(p, xs.*(1:n).*dx, ys.*pn./dx, label=string(i); kwargs...)
    end
    display(p)
    return p
end

plot_densities( ps, 11:16 )
plot_densities( ps, 11:16, factor=0.95 )
plot_densities( ps, 11:16, factor=0.95, yscale=:log10 )

p = plot_densities( ps, 16:16 )
p = plot_densities( ps, 16:16, xscale=0.5, p=p, yscale=0.25 )

M = calculate2( 0.9, 2 )

i = 11
factor=1.0
p = plot( size=[1000,800] )
xs=1.0
ys=1.0
kwargs = ()
pn = ps[i]
j = (f -> f == nothing ? length(pn) : f-1 )(findfirst( pn .== 0 ))
pn = pn[1:j]
n = length(pn)
dx = 1/((2 * factor) ^ i )
plot!(p, xs.*(1:n).*dx, ys.*pn./dx, label=string(i); kwargs...)
display(p)







