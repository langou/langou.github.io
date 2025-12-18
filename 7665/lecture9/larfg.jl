    function larfg!(x)
        """Computes the Householder transformation for input vector x"""

        σ = dot(x[2:end],x[2:end])

        if σ == 0
            τ = 0
            return τ
        end

        β = sqrt(x[1]^2 + σ)

        if x[1] > 0
            x[1] = x[1] + β
            x₁ = -β
        else
            x[1] = x[1] - β 
            x₁ = +β
        end

        x[2:end] = x[2:end] / x[1]

        τ = 2.0 / (1.0 + σ / x[1]^2)

        x[1] = x₁

        return τ

    end
