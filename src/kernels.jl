# X.F. Yang's Double Cosine kernel (2013)
function doublecosine(r, h)
    # Kernel constant
    k = 2
    if abs(r) < 2h
        # Normalizing factor
        σ = 1/(6k * h)

        # Checking sign
        sign = r == 0 ? 1 : r/abs(r)

        # Kernel computation
        W = (4cos(π/(2h) * r) + cos(π/h * r) + 3) * σ

        # Derivatives computation
        ∂x_W = (4π/(2h) * sin(π/(2h) * r) + π/h * sin(π/h * r)) * -sign * σ

        return (W,∂x_W)
    else
        return (0,0)
    end
    res
end

# J.J. Monaghan's Cubic Spline kernel (1992)
function cubicspline(r, h)
    if abs(r) < 2h
        # Reduced radius
        q = r/h

        # Normalizing factor
        σ = 2/(3h)

        # Checking sign
        sign = r == 0 ? 1 : r/abs(r)

        if abs(r) < h
            # Kernel computation
            W = (1 - 1.5q^2 + 0.75q^3) * σ

            # Derivatives computation
            ∂x_W = (2q - 1.5q^2) * sign/h^2

        else #if h <abs(r) < 2h
            # Kernel computation
            W = (2-q)^3/4 * σ

            # Derivatives computation
            ∂x_W = (q-2)^2/2 * sign/h^2
        end
        return (W,∂x_W)
    else
        return (0,0)
    end
    res
end
