
"""
    shadow_cylinder(pos::AbstractVector, sun::AbstractVector, R::Number)

Shadow function using cylindrical shadow model.

### Arguments 
- `occ` : observer to occulted body 
- `src` : observer to source body
- `R`   : occulted body radius 
"""
@fastmath function shadow_cylinder(occ::AbstractVector, src::AbstractVector, R::Number)   
    on = norm(occ)
    s = unitvec(src)

    # Compute the angle between the body-sun vector and the body-sun vector
    θ = acos( dot(-occ, s) / on )

    # The satellite is in eclipse if this angle is less than the 
    # planet's angular radius
    return if θ ≥ asin(R / on) 
        1.0
    else 
        0.0 
    end
end


"""
    shadow_conic(occ::AbstractVector, src::AbstractVector, R::Number, Rs::Number)

Shadow function using conical shadow model.

### Arguments 
- `occ` : observer to occulted body 
- `src` : observer to source body
- `R`   : occulted body radius 
- `Rs`  : source body radius
"""
@fastmath function shadow_conic(
    occ::AbstractArray, src::AbstractArray, R::Number, Rs::Number
)   
    on = norm(occ)
    sn = norm(src)

    a, b, c = asin(Rs/sn), asin(R/on), acos( dot(occ, src)/(on*sn) )
    a², b², c² = a^2, b^2, c^2

    if c ≥ (b + a)
        return 1.0
    elseif c < (b - a)
        return 0.0
    elseif c < (a - b)
        return 1.0 - b²/a²
    else
        x = (c² + a² - b²) / (2.0 * c)
        y = sqrt(a² - x^2)
        A = a² * acos(x / a) + b² * acos((c - x) / b) - c * y
        return 1.0 - A / (π * a²)
    end
end