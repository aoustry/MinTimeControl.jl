using Gurobi
using Random
using JuMP

include("system.jl")

########################################### System definition #####################################################
nx, nu = 2,2;
xm,um,xM,uM = -1.0,-1.0,1.0,1.0
Tmax = 10.0;
x0 = [0.0,-0.7];
xT = [0,0.7]
sys = zermelo_boat(nx,nu,0.5,xM*ones(2),uM*ones(2),xm*ones(2),um*ones(2));
##################################################################################################################
include("basis.jl");
include("solver.jl");
include("plots.jl");

trajec,tmax = heuristic_trajectory(sys,x0,xT,Tmax)


model = Model(Gurobi.Optimizer);
@variable(model, θ[1:N]);
constraints = [];
add_hamiltonian_constraints_on_trajectory(model, θ, sys, trajec,1,constraints);
add_hamiltonian_constraints_random(model, θ, sys,tmax,20000,constraints);
add_terminal_constraints(model, θ,xT,tmax,100,constraints);
add_terminal_sdp_constraints_random(model, θ,xT,tmax,100,constraints);
add_terminal_gradient_constraints(model, θ,sys,xT,tmax,100,constraints);

set_objective(model,θ,x0,1e-6);
optimize!(model);

λ = [value(θ[i]) for i in 1:N];
tr,time_feas = vmin_trajectory(sys,x0,xT,2.0,λ)

name_system = "zermelo"



plot_∇v_terminal(sys,λ,xT,tmax,name_system)
plot_v_terminal(sys,λ,xT,tmax,name_system)
plot_v(sys,λ,tmax,name_system);
plot_traj(tr,[2,3],name_system)
plot_v_traj(sys,λ,tr,name_system)


println("Time feasible control = ",time_feas)
println("Value v(0,x0) = ",v(vcat(0.0,x0),λ));

#ToBeDeleted : 
# Cycle
#Using findall(c->c==constraint[10],constraints)
#Using deleteat!(constraints, idx)