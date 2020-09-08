using LinearAlgebra

eye(N) = Matrix(1.0I,N,N)


##

# There are several root-finding routines :
# (1) python scipy
# (2) Julia NLsolve
# (3) Minpack (written in C, wrapped by Julia)

# WARNING the choice of root-finding routine has major impact
# on speed and reliability of the integrator

using NLsolve # reference: http://fourier.eng.hmc.edu/e176/lectures/NM/node21.html
# using MINPACK

function find_root(residue::Function, jacobian::Function, x0)
    sol = NLsolve.nlsolve(residue, jacobian, x0, method = :newton)
    return sol.zero
end

#using MINPACK

#function find_root(residue::Function, jacobian::Function, x0)
#    sol = NLsolve.nlsolve(residue, jacobian, x0, method = :newton)
#    return sol.zero
#end

#using PyCall
#sp = pyimport("scipy")

##

mutable struct Integrator
    H::Function
    ∇H∇p::Function          # must be kept for ALL types of integrators
    ∇H∇q::Function          # must be kept for ALL types of integrators
    ∇H∇p∇p::Function        # can be removed for certain types of integrators
    ∇H∇p∇q::Function        # must be kept for ALL types of integrators
    ∇H∇q∇q::Function        # can be removed for certain types of integrators
    p::Vector{Float64}
    q::Vector{Float64}
    h::Float64
end

@inline new_Integrator( IG::Integrator,
                        p_new::Vector,
                        q_new::Vector ) = Integrator( IG.H,
                                                      IG.∇H∇p,
                                                      IG.∇H∇q,
                                                      IG.∇H∇p∇p,
                                                      IG.∇H∇p∇q,
                                                      IG.∇H∇q∇q,
                                                      p_new, q_new,
                                                      IG.h  )


## -------------- INTEGRATORS ------------------
# REFERENCE  https://github.com/algorithmx/JuliaBasicExamples/blob/master/week2.pdf

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
    p1_sol = find_root(fp!, jp!, IG.p .+ (random_amp*IG.h).*rand(dim))
    q1     = IG.q .+ IG.h .* IG.∇H∇p(p1_sol, IG.q)
    return new_Integrator( IG, p1_sol, q1 )
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
    q1_sol = find_root(fq!, jq!, IG.q .+ (random_amp*IG.h).*rand(dim))
    p1     = IG.p .- IG.h .* IG.∇H∇q(IG.p, q1_sol)
    return new_Integrator( IG, p1, q1_sol )
end


# "the shorter, the faster"

@inline function Symplectic_Euler_V1_pq_separable(IG::Integrator)
    p1_sol = IG.p .- IG.h .* IG.∇H∇q(0,IG.q)
    return new_Integrator( IG, p1_sol, IG.q .+ IG.h .* IG.∇H∇p(p1_sol,0) )
end


@inline function Symplectic_Euler_V2_pq_separable(IG::Integrator)
    q1_sol = IG.q .+ IG.h .* IG.∇H∇p(IG.p,   0)
    return new_Integrator( IG, IG.p .- IG.h .* IG.∇H∇q(0,q1_sol), q1_sol )
end
