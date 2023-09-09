function facdiv(n::Int, m::Int)
    res = 1
    tmp = n
    while tmp > m
        res *= tmp
        tmp -= 1
    end
    return res
end

function normalization_factors(p::Int, q::Int, n::Int, m::Int)
    numerator = q == 0 ? 1 : 2
    denominator = m == 0 ? 1 : 2
    numerator *= 2p + 1
    denominator *= 2n + 1

    if n + m > p + q
        numerator *= facdiv(n + m, p + q)
    else
        denominator *= facdiv(p + q, n + m)
    end

    if p - q > n - m
        numerator *= facdiv(p - q, n - m)
    else
        denominator *= facdiv(n - m, p - q)
    end
    return sqrt(numerator / denominator)
end

@inline compute_η0(m) = (2m-1) * normalization_factors( m, m, m - 1, m - 1 )
@inline compute_η1(n, m) = (2*n-1)/(n-m) * normalization_factors( n, m, n - 1, m )
@inline compute_η2(n, m) = (n+m-1)/(n-m) * normalization_factors( n, m, n - 2, m )