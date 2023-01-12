using Gurobi
using Random
using JuMP

dt = 1e-3;
RADIUS_BALL = TARGET_TOLERANCE = 5*1e-2;
include("system.jl")
Random.seed!(0)
N_TERMINAL = 1000
N_SELECTED_CONSTRAINTS = 100 ;
MAX_NUMBER_OF_ATTEMPTS = 100000;
N_SELECTED_CONSTRAINTS_FINAL = 100 ;
MAX_NUMBER_OF_ATTEMPTS_FINAL = 10000;
N_RANDOM_INIT = 1000
CIRCLE_FINAL_SET = true;
TERMINAL_CONSTRAINT_EQ = false; 
μ = 1e-5;
STEP_ADD_TRAJ = 5
########################################### System definition #####################################################
#include("system_definition_zermelo.jl")
include("system_definition_toyboat.jl")
##################################################################################################################
include("basis.jl");
include("solver.jl");
include("plots.jl");


traj_heur,tmax = heuristic_trajectory(sys,x0,xT,Tmax)
best_traj,tmax,λ = dual_solving(sys,x0,xT,traj_heur,tmax,0.1,μ)
plot_∇v_terminal(sys,λ,xT,tmax)
plot_v_terminal(sys,λ,xT,tmax)
plot_traj(best_traj,[2,3])
plot_v_traj(sys,λ,best_traj)
plot_hmin_traj(sys,λ,best_traj)
#plot_hmin(sys,λ,tmax*1.1)
#plot_field(sys,λ,tmax*1.1)
