using Gurobi
using Random
using JuMP
using JLD
Random.seed!(0)
dt = 1e-3; #1e-3
RADIUS_BALL = TARGET_TOLERANCE = 5*1e-2; #1e-2
include("system.jl")
N_TERMINAL = 1000
N_SELECTED_CONSTRAINTS = 100 ;
MAX_NUMBER_OF_ATTEMPTS = 100000;
N_SELECTED_CONSTRAINTS_FINAL = 100 ;
MAX_NUMBER_OF_ATTEMPTS_FINAL = 10000;
N_RANDOM_INIT = 1000
CIRCLE_FINAL_SET = true;
TERMINAL_CONSTRAINT_EQ = false; 
μ = 1e-5;  #1e-5
ϵ = 1e-3;   #0.05
STEP_ADD_TRAJ = 1 ; #5
OBJECTIVE_WITH_INTEGRAL = false; 
simplex = true;
########################################### System definition #####################################################
#system_file_name = "system_definition_zermelo.jl";
#system_file_name ="system_definition_toyboat.jl";
system_file_name ="system_definition_gen_brockett6.jl";
include(system_file_name);
##################################################################################################################
include("basis.jl");
include("solver.jl");
include("plots.jl");
##################################################################################################################
time = @elapsed begin
success,traj_heur,tmax = heuristic_trajectory(sys,x0,xT,Tmax)
if success==false
    tmax = Tmax
end
println("Success heuristic trajectory = ",success)
println("Tmax heuristic trajectory = ",tmax)
best_traj,tmax,λ, lb, iter, seq_control,seq_ub,seq_lb,seq_bounds = dual_solving(sys,x0,xT,traj_heur,tmax,ϵ,μ,simplex)
end

res =  Dict("time"=>time,"ub"=> tmax,"lb"=> lb,"degree"=> degree,"iterations"=> iter,"array_time_control"=>seq_control,"array_ub_siprho"=>seq_ub,"array_lb_siprho"=>seq_lb,"array_lb"=>seq_bounds)

array = split(system_file_name,"_")
name = array[length(array)]
save("logs/"*name*"_"*string(degree)*".jld", "data", res)

#plot_∇v_terminal(sys,λ,xT,tmax)
#plot_v_terminal(sys,λ,xT,tmax)
#plot_traj(best_traj,[2,3])
#plot_v_traj(sys,λ,best_traj)

