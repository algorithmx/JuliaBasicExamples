include("integrators.jl")

using BenchmarkTools

using Plots

using StatsBase

using JLD2


# automatic differentiation
# either to verify the gradients / hessians
# or to get them automatically
using ForwardDiff


## 

# particles in a box with Lennard-Jones interaction (3D)

# settings

N = 10^3   # number of particles
L = 1.0    # size of the bounding box

m = 0.2    # mass of the particle  
minv = 1.0/m
m2inv = 0.5/m
ε = 3.0    # strength parameter of LJ-6-12 potential
rm = 0.02  # r_min  parameter of LJ-6-12 potential
rm2 = rm^2
rm6 = rm^6

# auxiliary matrices

Minv = (1.0/m).*ones(Float64,3N)
M2inv = (0.5/m).*ones(Float64,3N)

zero3Nx3N = zeros(Float64, (3N,3N))

X = reshape(zeros(Float64,3,N).+[1,0,0],3N)
Y = reshape(zeros(Float64,3,N).+[0,1,0],3N)
Z = reshape(zeros(Float64,3,N).+[0,0,1],3N)

XX = reshape(zeros(Float64,3,N).+[400,0,0],3N)
YY = reshape(zeros(Float64,3,N).+[0,400,0],3N)
ZZ = reshape(zeros(Float64,3,N).+[0,0,400],3N)

mXX = reshape(zeros(Float64,3,N).+[-400,0,0],3N)
mYY = reshape(zeros(Float64,3,N).+[0,-400,0],3N)
mZZ = reshape(zeros(Float64,3,N).+[0,0,-400],3N)

XX2 = reshape(zeros(Float64,3,N).+[400^2,0,0],3N)
YY2 = reshape(zeros(Float64,3,N).+[0,400^2,0],3N)
ZZ2 = reshape(zeros(Float64,3,N).+[0,0,400^2],3N)


## ------ confining potential ------------

U_confine(q) = ( sum(X.*exp.(mXX.*q) .+ X.*exp.((XX.*q).-XX)) + 
                 sum(Y.*exp.(mYY.*q) .+ Y.*exp.((YY.*q).-YY)) +
                 sum(Z.*exp.(mZZ.*q) .+ Z.*exp.((ZZ.*q).-ZZ)) )

dU_confine(q)  = ( (mXX.*exp.(mXX.*q) .+ XX.*exp.((XX.*q).-XX)) .+ 
                   (mYY.*exp.(mYY.*q) .+ YY.*exp.((YY.*q).-YY)) .+
                   (mZZ.*exp.(mZZ.*q) .+ ZZ.*exp.((ZZ.*q).-ZZ)) )

ddU_confine(q) = diagm( 0=>(XX2.*(exp.(mXX.*q) .+ exp.((XX.*q).-XX)) .+ 
                            YY2.*(exp.(mYY.*q) .+ exp.((YY.*q).-YY)) .+
                            ZZ2.*(exp.(mZZ.*q) .+ exp.((ZZ.*q).-ZZ))) )

# test
ddU_confine1(q)  = ForwardDiff.hessian(q->U_confine(q), q)
dU_confine1(q)   = ForwardDiff.gradient(q->U_confine(q), q)
#@assert all([(q=rand(3N); dU_confine(q) ≈ dU_confine1(q)) for i=1:3])
#@assert all([(q=rand(3N); ddU_confine(q) ≈ ddU_confine1(q)) for i=1:3])


## ------ LJ-6-12 ------------

# for a pair of particles

to6(dR) = (rm2/dot(dR,dR))^3
Rn(dR) = dR./dot(dR,dR)

V_LJ_6_12(rm6_r6) = ε*((rm6_r6-1)^2-1)
dV_LJ_6_12_X(rm6_r6, rn) = (12ε*rm6_r6*(1-rm6_r6)) .* rn
dV_LJ_6_12(dr) = dV_LJ_6_12_X(to6(dr), Rn(dr))

# test
dV_LJ_6_12_a(r1,r2) = ForwardDiff.gradient(r1->V_LJ_6_12(to6(r1.-r2)), r1)
r1 = 1.0.+0.1rand(3)
r2 = 1.1.+0.1rand(3)
dV1 = dV_LJ_6_12( r1.-r2 )
dV2 = dV_LJ_6_12_a( r1, r2 )
@assert dV1 ≈ dV2

# for all particles

U_LJ(q) = sum( V_LJ_6_12(to6(q[3i-2:3i].-q[3j-2:3j])) for i=1:N for j=1:i-1 )

function dU_LJ(q)
    dU = zeros(Float64,length(q))
    for i=1:N
        qi = q[3i-2:3i]
        dU[3i-2:3i] .= sum( dV_LJ_6_12(qi.-q[3j-2:3j]) for j=1:N if i!=j )
    end
    dU
end

ddU_LJ(q) = ForwardDiff.hessian(q->U_LJ(q), q)  # I'm lazy

#test
dU_LJ_a(q) = ForwardDiff.gradient(q->U_LJ(q), q)
#@assert all([(q = 2rand(3N); dU_LJ(q)≈dU_LJ_a(q)) for i = 1:3])

##

# TODO speed up !!!
                
U(q)   = U_confine(q) + U_LJ(q)
dU(q)  = dU_confine(q) .+ dU_LJ(q)
ddU(q) = ddU_confine(q) .+ ddU_LJ(q)

T(p)   = m2inv.*(p.*p)
dT(p)  = minv.*p
ddT(p) = Minv

# canonical momentum p and coordinate q are separable
H(p,q)      = T(p) + U(q)
∇H∇p(p,q)   = dT(p)
∇H∇q(p,q)   = dU(q)
∇H∇q∇q(p,q) = ddU(q)
∇H∇p∇p(p,q) = ddT(p)
∇H∇p∇q(p,q) = zero3Nx3N

## ------------------ MD, main ---------------------

# computation takes long time

#=

# N = 10^3
Q0 = vcat([-0.05.+0.005rand(3).+[0.1i,0.1j,0.1k] for i=1:10 for j=1:10 for k=1:10]...)
P0 = rand(3N)
IG = Integrator( H, ∇H∇p, ∇H∇q, ∇H∇p∇p, ∇H∇p∇q, ∇H∇q∇q, P0, Q0, 0.0001 ) ;

##

Nsteps = 8000

traj = [IG,]

for i = 1:Nsteps
    if i%100==0
        println("Progress ", round(i/Nsteps,digits=3),"%")
    end
    IG1 = Symplectic_Euler_V1(traj[end])
    push!(traj, IG1)
end

@save "traj.jld2" traj

=#

## --------------------------------------------------

## load results (computation takes long time)

@load "traj.jld2"
Q = hcat([t.q for t ∈ traj]...) ;
P = hcat([t.p for t ∈ traj]...) ;
Qxyz = reshape(Q,3,N,:) ;
Pxyz = reshape(P,3,N,:) ;

## ----------- pair-distribution function -----------

# TODO

Nsteps = 8000

dist(Q, r0) = vec(sqrt.(sum((Q.-r0).^2,dims=1)))

findr0(Q) = Q[:,findfirst(dist(Q,[0.5,0.5,0.5]).<0.1)]

in_shell_0(d, R, dR) = findall((d.<(R+dR)).&(d.>=R))
in_shell(Q, R, dR) = in_shell_0(dist(Q,findr0(Q)), R, dR)

PDF(Q, dr, rmax) = [(1.0/(4π*1000*((R+0.5dr)^2)))*length(in_shell(Q,R,dr)) for R=(rm-1e-10):dr:rmax]

pdf = hcat([PDF(Qxyz[:,:,k], 0.01, 1) for k=2000:20:Nsteps]...)

##

pdf_mean = vec(mean(pdf,dims=2))

plot( collect(rm-1e-10:0.01:1)[1:40], (pdf_mean.-1)[1:40] )



## ----------- position animations ---------------

anim = @animate for i = 3000:20:Nsteps
    scatter(Qxyz[1,:,i],Qxyz[2,:,i], 
            aspect_ratio=1, xlims=(-0.05,1.05), ylims=(-0.05,1.05))
end

gif(anim,  "XY_plane.gif", fps = 15)


anim = @animate for i = 3000:20:Nsteps
    scatter(Qxyz[2,:,i],Qxyz[3,:,i], 
            aspect_ratio=1, xlims=(-0.05,1.05), ylims=(-0.05,1.05))
end

gif(anim,  "YZ_plane.gif", fps = 15)


anim = @animate for i = 3000:20:Nsteps
    scatter(Qxyz[1,:,i],Qxyz[3,:,i], 
            aspect_ratio=1, xlims=(-0.05,1.05), ylims=(-0.05,1.05))
end

gif(anim,  "XZ_plane.gif", fps = 15)


## ----------- momentum statistics --------------

# TODO ???

Pn = [dot(Pxyz[:,i,k],Pxyz[:,i,k]) for i=1:N for k=3000:20:Nsteps]
hist_pdf = normalize(fit(Histogram, Pn, 0:0.1:3.0))
plt = plot(hist_pdf.edges[1][1:end-1], hist_pdf.weights)


## show equipartition

PnX = [abs(Pxyz[1,i,k])^2 for i=1:N for k=3000:20:Nsteps]
PnY = [abs(Pxyz[2,i,k])^2 for i=1:N for k=3000:20:Nsteps]
PnZ = [abs(Pxyz[3,i,k])^2 for i=1:N for k=3000:20:Nsteps]
hist_pdf_X = normalize(fit(Histogram, PnX, 0:0.1:3.0))
hist_pdf_Y = normalize(fit(Histogram, PnY, 0:0.1:3.0))
hist_pdf_Z = normalize(fit(Histogram, PnZ, 0:0.1:3.0))
plt = plot(hist_pdf_X.edges[1][1:end-1], hist_pdf_X.weights)
plt = plot!(hist_pdf_Y.edges[1][1:end-1], hist_pdf_Y.weights)
plt = plot!(hist_pdf_Z.edges[1][1:end-1], hist_pdf_Z.weights)
plt


## -------------------------------------------------------------------------------------
# accelerated for the separable Hamiltonian H = T(p) + U(q)

#=

# N = 10^3
Q0 = vcat([-0.05.+0.005rand(3).+[0.1i,0.1j,0.1k] for i=1:10 for j=1:10 for k=1:10]...)
P0 = rand(3N)
IG = Integrator( H, ∇H∇p, ∇H∇q, ∇H∇p∇p, ∇H∇p∇q, ∇H∇q∇q, P0, Q0, 0.0001 ) ;

Nsteps = 8000

traj = [IG,]

for i = 1:Nsteps
    if i%100==0
        println("Progress ", round(i/Nsteps,digits=3),"%")
    end
    IG1 = Symplectic_Euler_V1_pq_separable(traj[end])
    push!(traj, IG1)
end

@save "traj_1.jld2" traj

=#