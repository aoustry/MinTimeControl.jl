using Gurobi
using Random
using JuMP

include("system.jl")
Random.seed!(0)
N_TERMINAL = 1000
N_SELECTED_CONSTRAINTS = 100 ;
N_RANDOM_INIT = 1000
########################################### System definition #####################################################
#include("system_definition_zermelo.jl")
include("system_definition_toyboat.jl")
##################################################################################################################
include("basis.jl");
include("solver.jl");
include("plots.jl");


function dual_solving(sys,x0,xT,traj_heur,tmax,ϵ)
    best_traj = traj_heur
    model = Model(Gurobi.Optimizer);
    set_optimizer_attribute(model, "Method", 1)
    @variable(model, θ[1:N]);
    constraints = [];
    add_hamiltonian_constraints_on_trajectory(model, θ, sys, traj_heur,constraints);
    add_hamiltonian_constraints_random(model, θ, sys,tmax,N_RANDOM_INIT,constraints);
    add_terminal_constraints(model, θ,xT,tmax,N_TERMINAL,constraints);
    add_terminal_gradient_constraints(model, θ,sys,xT,tmax,N_TERMINAL,constraints);
    add_terminal_sdp_constraints(model, θ,xT,tmax,200,constraints);

    set_objective(model,θ,x0,1e-5);
    λ = 0.0;
    max_violation = -Inf;
    while max_violation < -ϵ
        optimize!(model);
        λ = [value(θ[i]) for i in 1:N];
        traj,tmax_new = vmin_trajectory(sys,x0,xT,tmax,λ)
        if tmax_new < tmax
            tmax = tmax_new
            best_traj = traj
        end
        println("Time best feasible control = ",tmax);
        println("Value v(0,x0) = ",v(vcat(0.0,x0),λ));
        add_hamiltonian_constraints_on_trajectory(model, θ, sys, traj,constraints,λ);
        max_violation = add_selected_cuts(model, θ, sys,tmax,N_SELECTED_CONSTRAINTS,constraints,λ)
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
    return best_traj, tmax,λ
end


traj_heur,tmax = heuristic_trajectory(sys,x0,xT,Tmax)
best_traj,tmax,λ = dual_solving(sys,x0,xT,traj_heur,tmax,0.1)
plot_∇v_terminal(sys,λ,xT,tmax)
plot_v_terminal(sys,λ,xT,tmax)
plot_traj(best_traj,[2,3])
plot_v_traj(sys,λ,best_traj)
plot_hmin(sys,λ,tmax*1.1)
plot_field(sys,λ,tmax*1.1)



#ToBeDeleted : 
# Keep track of multiple trajectories
# plot Hmin
#Using findall(c->c==constraint[10],constraints)
#Using deleteat!(constraints, idx)