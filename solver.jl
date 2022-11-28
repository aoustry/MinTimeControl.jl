using Gurobi
using JuMP
using Random


nx, nu = 2,2;
Tmax = 10.0;
x0 = [0.0,-0.5];
xT = [0,0.5]
sys = zermelo_boat(nx,nu,0.8,ones(2),ones(2),-ones(2),-ones(2));

include("system.jl")
include("basis.jl")

function set_objective(model,θ,x0::Vector{Float64})
    a = ϕ(vcat(0,x0));
    @objective(model, Max, sum(a[i]*θ[i] for i in 1:N));
end

function add_hamiltonian_constraint(model,θ, sys::system,t::Float64,x::Vector{Float64},u::Vector{Float64})
    a =  ∇ϕ(vcat(t,x)) *vcat(1,f(sys,t,x,u));
    @constraint(model,sum(a[i]*θ[i] for i in 1:N)>=-1);
end

function add_terminal_constraint(model, θ, t::Float64, xT::Vector{Float64})
    a = ϕ(vcat(t,xT));
    @constraint(model,sum(a[i]*θ[i] for i in 1:N)<=0);
end

function add_terminal_sdp_constraint(model, θ, t::Float64, xT::Vector{Float64},g::Vector{Float64})
    a = g∇²ϕg(vcat(t,xT),g);
    @constraint(model,sum(a[i]*θ[i] for i in 1:N)=>0);
end

trajec,tmax = heuristic_trajectory(sys,x0,xT,Tmax)
model = Model(Gurobi.Optimizer);
@variable(model, θ[1:N]);

for vector in trajec
    add_hamiltonian_constraint(model,θ, sys,vector[1],vector[1:nx+1],vector[nx+2:nx+nu+1]);
end

timing = @elapsed begin
    for k in 1:1000
        t = 0.1
        x = [Random.rand(),Random.rand()]
        u = [1.0,0.0]
        add_hamiltonian_constraint(model,θ, sys,t,x,u);
    end
end

#ToBeDeleted:#Faire un fichier auxiliaire pour les fonctions types

#delete a constraint
#delete(model,constraint_by_name(model,"name"))