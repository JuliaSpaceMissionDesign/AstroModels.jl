
# Compute parametrized post newtonian correction 
# bodies[1] is the target body, bodies[2:end] the considered one for the correction 
# pcache, vcache and acache are caches of the barycentric position, velocity and accelerations of bodies 
# 
# Implemented according to Mayer, Formulation for Observed and Computed Values of Deep Space 
# Network Data Types for Navigation. Equation (4-26) without the point mass term

function compute_relativistic_ppn(pcache, vcache, acache, bodies, β, γ, c)
    # TODO: some computations here can be cached (Ul, Uk, ṡj and rji/rli/rjk computations)
    c2 = c*c
    a = SVector{3}(0, 0, 0) 
    Vi = vcache[1]
    ṡi = norm(Vi)

    for j in eachindex(bodies) if j != 1

        μj = bodies[j]
        Rji = pcache[j] - pcache[1]
        rji = norm(Rji)

        Vj = vcache[j]
        Aj = acache[j]

        Rji_Vj = dot(Rji, Vj)
        ṡj = norm(Vj)

        Ul = 0 
        for l in eachindex(bodies) if l != 1 
            rli = norm(pcache[l] - pcache[1]) 
            Ul += bodies[l].μ/rli
        end end

        Uk = 0 
        for k in eachindex(bodies) if k != l 
            rjk = norm(pcache[j] - pcache[k]) 
            Uk += bodies[k].μ/rjk
        end end

        F1 = -2*(β+γ)/c2*Ul -(2β-1)/c2*Uk + γ*ṡi^2/c2 + (1+γ)*ṡj^2/c2 - 2*(1+γ)/c2*dot(Vi, Vj)
        F1 += 3/(2*c2)*( Rji_Vj/rji )^2 + 1/(2*c2)*dot(Rji, Aj)
        F2 = - (2+2γ)*dot(Rji, Vi) + (1+2γ)*Rji_Vj

        a += μj/rji^3 * ( F1 * Rji + F2 * (Vi - Vj) ) + (3+4γ)/(2*c2)*μj/rji * Aj

    end end 
    return a
end

struct RelativisticPPNCorrections{T} <: AbstractGravityModel{T}
    bodies::Vector{PointMass{T}}
    β::T 
    γ::T 
    c::T

    # caches
    pos::Vector{SVector{3, T}}
    vel::Vector{SVector{3, T}}
    acc::Vector{SVector{3, T}}
end

function RelativisticPPNCorrections(bodies::Vector{PointMass{T}}, β, γ, c)
    nbodies = length(bodies)
    pc = [@SVector zeros(T, 3) for _ in 1:nbodies]
    vc = similar(pc)
    ac = similar(pc)
    return RelativisticPPNCorrections{T}(bodies, β, γ, c, pc, vc, ac)
end

function compute_relativistic_ppn(rc::RelativisticPPNCorrections{T}, frames::G, epoch) where {T, G}
    # precompute terms 
    @inbounds for i in eachindex(bodies) 
        b = bodies[i]
        pva = vector12(frames, 0, b.id, 1, epoch)
        rc.pos[i] = @views pva[1:3]
        rc.vel[i] = @views pva[4:6]
        rc.acc[i] = @views pva[7:9]
    end
    return compute_relativistic_ppn(rc.pos, rc.vel, rc.acc, rc.bodies, rc.β, rc.γ, rc.c)
end