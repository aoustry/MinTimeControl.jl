using Gurobi
using Random
using JuMP
using JLD
dt = 1e-3; #1e-3
RADIUS_BALL = TARGET_TOLERANCE = 5*1e-2; #1e-2
include("system.jl")
N_TERMINAL = 1000
N_SELECTED_CONSTRAINTS = 100 ;
MAX_NUMBER_OF_ATTEMPTS = 500000;
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
system_file_name ="system_definition_toyboat.jl";
#system_file_name ="system_definition_gen_brockett6.jl";
include(system_file_name);
##################################################################################################################
include("basis.jl");
include("solver.jl");
include("plots.jl");
##################################################################################################################

n = 1+nx;
@polyvar xvar[1:n];

function main(degree,certification)

    
    global basis = monomials(xvar, 0:degree)
    global N = length(basis)
    global symb∇ = [(differentiate(p, xvar)) for p in basis]


    time_heur = @elapsed begin
    success_heur,traj_heur,tmax = heuristic_trajectory(sys,x0,xT,Tmax)
    if success_heur==false
        tmax = Tmax
    end
    end
    println("Success heuristic trajectory = ",success_heur)
    println("Tmax heuristic trajectory = ",tmax)
    Random.seed!(0)
    best_traj,tmax,λ, lb, logs = dual_solving(sys,x0,xT,traj_heur,tmax,ϵ,μ,simplex,certification)

    logs["time_heur"] = time_heur
    logs["ub"] = tmax
    logs["lb"] = lb
    #Saving result
    array = split(system_file_name,"_")
    name = array[length(array)]
    save("logs/"*name*"_"*string(degree)*"_"*string(certification)*".jld", "data", logs)
    #plot_control_zermelo(sys,traj_heur,best_traj)
    #plot_control(sys,traj_heur,best_traj)
    #plot_water_flow(sys,traj_heur,best_traj)
    #plot_control_regatta(sys,traj_heur,best_traj)
end

certif = false
main(6,certif);

#for i in 2:8
#    main(i,certif);
#end
 
 

