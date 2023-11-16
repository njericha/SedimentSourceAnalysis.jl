# This file is to test the idea of simutaniously decomposing many multivariate distributions
# into a linear combination of a few multivaritate distributions via tensor decompositions.
#
# We do not assume the distributions are product distributions.

using Random

"""
Returns a sample [x1, x2] from a distribution B1 where
- x1 is standard normal
- x2 = x1 + z/10 and z is standard normal
"""
function B1()
    x1 = randn()
    x2 = x1 + randn()/10
    return [x1, x2]
end

using Plots

# Plot some samples
samples = [B1() for _ in 1:25]

x1 = [x for (x, _) in samples]
x2 = [x for (_, x) in samples]

scatter(x1, x2)

"""
Returns a sample [x1, x2] from a distribution B1 where
- x1 is standard normal
- x2 = -x1^2 + 1 + z/10 and z is standard normal
"""
function B2()
    x1 = randn()
    x2 = -x1^2 + 1 + randn()/10
    return [x1, x2]
end

# Plot some more samples
samples = [B2() for _ in 1:25]

x1 = [x for (x, _) in samples]
x2 = [x for (_, x) in samples]

scatter(x1, x2)

# Combine into a sample that is a sample from B1 with prob p and B2 with prob 1-p
"""
    B(p)

Samples B1 with probability p, and B2 with probability (1-p). p must be positive.
"""
function B(p)
    roll = rand()
    if roll <= p
        return B1()
    else
        return B2()
    end
end

# Plot some more samples
samples = [B(0.5) for _ in 1:100]

x1 = [x for (x, _) in samples]
x2 = [x for (_, x) in samples]

scatter(x1, x2)


# Perform 2dKDE
