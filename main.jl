using Gurobi
using Random
using JuMP

include("system.jl")
Random.seed!(0)
N_TERMINAL = 100
N_RANDOM_INIT = 1000
########################################### System definition #####################################################
param = 0.7
name_system = "zermelo"*string(param)
nx, nu = 2,2;
xm,um,xM,uM = -1.0,-1.0,1.0,1.0
Tmax = 10.0;
x0 = [0,-0.7];
xT = [0,0.7]
sys = zermelo_boat(nx,nu,param,xM*ones(2),uM*ones(2),xm*ones(2),um*ones(2));
##################################################################################################################
include("basis.jl");
include("solver.jl");
include("plots.jl");

traj_heur,tmax = heuristic_trajectory(sys,x0,xT,Tmax)
best_traj = traj_heur

model = Model(Gurobi.Optimizer);
set_optimizer_attribute(model, "Method", 1)
@variable(model, θ[1:N]);
constraints = [];
add_hamiltonian_constraints_on_trajectory(model, θ, sys, traj_heur,1,constraints);
add_hamiltonian_constraints_random(model, θ, sys,tmax,N_RANDOM_INIT,constraints);
add_terminal_constraints(model, θ,xT,tmax,N_TERMINAL,constraints);
add_terminal_gradient_constraints(model, θ,sys,xT,tmax,N_TERMINAL,constraints);
#add_terminal_sdp_constraints(model, θ,xT,tmax,200,constraints);

set_objective_with_integral(model,θ,x0,1e-3,tmax);
λ = 0.0;
max_violation = -Inf;
while max_violation < -0.1
    optimize!(model);
    λ = [value(θ[i]) for i in 1:N];
    traj,tmax_new = vmin_trajectory(sys,x0,xT,tmax,λ)
    if tmax_new < tmax
        tmax = tmax_new
        best_traj = traj
    end
    println("Time best feasible control = ",tmax);
    println("Value v(0,x0) = ",v(vcat(0.0,x0),λ));
    add_hamiltonian_constraints_on_trajectory(model, θ, sys, traj,4,constraints,λ);
    max_violation = add_selected_cuts(model, θ, sys,tmax,100,constraints,λ)
end

optimize!(model);
λ = [value(θ[i]) for i in 1:N];
traj,tmax_new = vmin_trajectory(sys,x0,xT,tmax,λ)
if tmax_new < tmax
    tmax = tmax_new
    best_traj = traj
end

println("Time feasible control = ",tmax)
println("Value v(0,x0) = ",v(vcat(0.0,x0),λ));

plot_∇v_terminal(sys,λ,xT,tmax,name_system)
plot_v_terminal(sys,λ,xT,tmax,name_system)
plot_traj(best_traj,[2,3],name_system)
plot_v_traj(sys,λ,best_traj,name_system)


#ToBeDeleted : 
# plot Hmin
#Using findall(c->c==constraint[10],constraints)
#Using deleteat!(constraints, idx)