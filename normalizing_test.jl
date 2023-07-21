using BenchmarkingTools

# Test which version is faster

function repeat_normalizing(D, v)
    I, _, K = size(D)
    scales = repeat(v', I, 1, K)
    return D ./ scales
end

function loop_normalizing(D, v)
    D_copy = copy(D)
    for (slice, x) in zip(eachslice(D_copy,dims=(1,3)), v)
        slice ./= x
    end
    return D_copy
end

@benchmark repeat_normalizing(D,v) setup=(D=randn(50,50,50),v=randn(50))

@benchmark loop_normalizing(D,v) setup=(D=randn(50,50,50),v=randn(50))

# Looks like the loop is faster!
