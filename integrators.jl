using LinearAlgebra
using NLsolve
using Plots

eye(N) = Matrix(1.0I,N,N)

##

mutable struct Integrator
    H::Function
    ∇H∇p::Function
    ∇H∇q::Function
    ∇H∇p∇p::Function
    ∇H∇p∇q::Function
    ∇H∇q∇q::Function
    p::Vector{Float64}
    q::Vector{Float64}
    h::Float64
end


##

## order 1

function Symplectic_Euler_V1(
    IG::Integrator;
    random_amp = 0
    )::Integrator
    dim = length(IG.p)
    function fp!(Fp, p1)
        Fp .= (p1 .- IG.p .+ IG.h.*IG.∇H∇q(p1, IG.q))
    end
    function jp!(Jp, p1)
        Jp .= (eye(dim) .+ IG.h.*IG.∇H∇p∇q(p1, IG.q)')
    end
    sol = nlsolve(fp!, jp!, IG.p .+ (random_amp*IG.h).*rand(dim))
    p1_sol = sol.zero
    q1 = IG.q .+ IG.h .* IG.∇H∇p(p1_sol, IG.q)
    return Integrator(  IG.H, IG.∇H∇p, IG.∇H∇q, IG.∇H∇p∇p, IG.∇H∇p∇q, IG.∇H∇q∇q, 
                        p1_sol, q1, IG.h  )
end


function Symplectic_Euler_V2(
    IG::Integrator;
    random_amp = 0
    )::Integrator
    dim = length(IG.p)
    function fq!(Fq, q1)
        Fq .= (q1 .- IG.q .- IG.h.*IG.∇H∇p(IG.p, q1))
    end
    function jq!(Jq, q1)
        Jq .= (eye(dim) .- IG.h.*IG.∇H∇p∇q(IG.p, q1))
    end
    sol = nlsolve(fq!, jq!, IG.q .+ (random_amp*IG.h).*rand(dim))
    q1_sol = sol.zero
    p1 = IG.p .- IG.h .* IG.∇H∇q(IG.p, q1_sol)
    return Integrator(  IG.H, IG.∇H∇p, IG.∇H∇q, IG.∇H∇p∇p, IG.∇H∇p∇q, IG.∇H∇q∇q, 
                        p1, q1_sol, IG.h  )
end


## -------------------------------------------------
## tests
# one particle moving in a 3d parabolic, isotropic confining potential

I3 = eye(3)

U(q,M) = 0.5 * (q'*M*q)
dU(q,M) = (M*q)
ddU(q,M) = M

T(p,invM) = 0.5 * (p'*invM*p)
dT(p,invM) = (invM*p)
ddT(p,invM) = invM

H(p,q)      = T(p, I3) + U(q, 5.0.*I3)
∇H∇p(p,q)   = dT(p, I3)
∇H∇q(p,q)   = dU(q, 5.0.*I3)
∇H∇q∇q(p,q) = ddU(q, 5.0.*I3)
∇H∇p∇p(p,q) = ddT(p, I3)
∇H∇p∇q(p,q) = 0.0.*I3

##

IG = Integrator( H, ∇H∇p, ∇H∇q, ∇H∇p∇p, ∇H∇p∇q, ∇H∇q∇q, 
                 [5.0,0.2,0.0], [0.0,1.0,2.0], 0.01 ) ;

##

traj1 = [IG,]
for i = 1:4000
    IG1 = Symplectic_Euler_V1(traj1[end])
    push!(traj1, IG1)
end

##

traj2 = [IG,]
for i = 1:4000
    IG2 = Symplectic_Euler_V2(traj2[end])
    push!(traj2, IG2)
end

##

plot_traj(t,view) = plot(map(x->x.q[view[1]],t), map(x->x.q[view[2]],t), xlims=(-4,4), ylims=(-4,4))
plot_traj!(t,view) = plot!(map(x->x.q[view[1]],t), map(x->x.q[view[2]],t), xlims=(-4,4), ylims=(-4,4))

plot_traj(traj1,[3,2])
plot_traj!(traj2,[3,2])
