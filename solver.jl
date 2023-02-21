using JuMP


############################################# Objective ####################################################
function set_objective(model,θ,x0::Vector{Float64},ρ::Float64)
    a = ϕ(vcat(0,x0));
    if ρ>0
        @objective(model, Max, sum(a[i]*θ[i] for i in 1:N) - ρ*sum(θ[i]*θ[i] for i in 1:N));
    else
        @objective(model, Max, sum(a[i]*θ[i] for i in 1:N));
    end
end

function set_objective_with_integral(model,θ,x0::Vector{Float64},ρ::Float64,tmax::Float64,param::Float64)
    a = ϕ(vcat(0,x0));
    P = 10000;
    moments = zeros(N);
    r = 0.05
    for aux in 1:P
        t = 0;
        x = (((2* r) .* [Random.rand() for j in 1:nx]) .- r .+ x0);
        moments = moments .+ ϕ(vcat(t,x))/P;
    end

    if ρ>0.0
        @objective(model, Max, sum(a[i]*θ[i] for i in 1:N) + param*sum(moments[i]*θ[i] for i in 1:N)  - ρ*sum(θ[i]*θ[i] for i in 1:N));
    else
        @objective(model, Max, sum(a[i]*θ[i] for i in 1:N) + param*sum(moments[i]*θ[i] for i in 1:N));
    end
end
############################################################################################################

########################################## Individual constraints ##########################################
function add_hamiltonian_constraint(model,θ, sys::system,t::Float64,x::Vector{Float64},u::Vector{Float64})
    a =  ∇ϕ(vcat(t,x)) *vcat(1,f(sys,t,x,u));
    return @constraint(model,sum(a[i]*θ[i] for i in 1:N)>=-1);
end

function add_terminal_constraint(model, θ, t::Float64, xT::Vector{Float64})
    a = ϕ(vcat(t,xT));
    if TERMINAL_CONSTRAINT_EQ
        return @constraint(model,sum(a[i]*θ[i] for i in 1:N)==0);
    else
        return @constraint(model,sum(a[i]*θ[i] for i in 1:N)<=0);
    end
end

function add_terminal_sdp_constraint(model, θ, t::Float64, xT::Vector{Float64},g::Vector{Float64})
    a = g∇²ϕg(vcat(t,xT),vcat(0,g));
    return @constraint(model,sum(a[i]*θ[i] for i in 1:N)>=0);
end

#function add_terminal_gradient_constraint(model, θ,sys::system, t::Float64, xT::Vector{Float64})
#    ∇ = ∇ϕ(vcat(t,xT));
#    mat = ∇*Diag(vcat(0,ones(sys.nx)))*∇' + 1e-6*I;
#    return @constraint(model,θ'*mat*θ<=1e-3);
#end

function add_terminal_gradient_constraint(model, θ,sys::system, t::Float64, xT::Vector{Float64})
    ∇ = ∇ϕ(vcat(t,xT));
    return @constraint(model,(∇'*θ)[2:1+sys.nx] .==0);
end

############################################################################################################


########################################## Multiple constraints ##########################################

function add_terminal_constraints(model, θ,xT::Vector{Float64},tmax::Float64,P::Int64,constraints::Vector{Any})
    step = tmax/P;
    for t in 0:step:tmax
        append!(constraints, [add_terminal_constraint(model,θ, t,xT)]);
    end
end

function add_terminal_sdp_constraints(model, θ,xT::Vector{Float64},tmax::Float64,P::Int64,constraints::Vector{Any})
    for t in tmax*0.5:tmax/(2*P):tmax
        for k in 1:20
            g = random_control(sys,t,xT);
            append!(constraints, [add_terminal_sdp_constraint(model,θ, t,xT,g)]); 
        end
    end
end

function add_terminal_gradient_constraints(model, θ,sys::system,xT::Vector{Float64},tmax::Float64,P::Int64,constraints::Vector{Any})
    for t in tmax*0.5:tmax/(2*P):tmax
        append!(constraints, [add_terminal_gradient_constraint(model,θ,sys,t,xT)];)
    end
end

function add_hamiltonian_constraints_on_trajectory(model, θ, sys::system, traj::Vector{Any},constraints::Vector{Any},λ::Vector{Float64})
    added = 0
    for aux in 1:STEP_ADD_TRAJ:length(traj)
            vector = traj[aux]; 
            t,x,u = vector[1],vector[2:nx+1],vector[nx+2:nx+nu+1];
            if is_feasible(sys,x) && H(sys,t,x,u,λ)<-1
                added+=1;
                append!(constraints, [add_hamiltonian_constraint(model,θ, sys,t,x,u)]);
            end
    end
    println("Added cuts along trajectory = ",added);
end


function add_hamiltonian_constraints_on_trajectory(model, θ, sys::system, traj::Vector{Any},constraints::Vector{Any})
    for aux in 1:length(traj)
            vector = traj[aux]; 
            x = vector[2:nx+1];
            if is_feasible(sys,x)
                append!(constraints, [add_hamiltonian_constraint(model,θ, sys,vector[1],x,vector[nx+2:nx+nu+1])]);
            end
    end
end

function add_hamiltonian_constraints_random(model, θ, sys::system,tmax::Float64,P::Int64,constraints::Vector{Any})
   for aux in 1:P
        t = Random.rand()*tmax;
        x = ((sys.xmax.-sys.xmin) .* [Random.rand() for j in 1:sys.nx]) .+ sys.xmin;
        if is_feasible(sys,x)
            u = random_control(sys,t,x);
            append!(constraints, [add_hamiltonian_constraint(model,θ, sys,t,x,u)]);
        end
   end
end

function random_vector_unit_ball(p::Int64)
    while true
        u = 2* [Random.rand() for j in 1:p] .- 1;
        if norm(u)<=1
        return u
        end
    end
end
############################################################################################################


########################################## Random search function ##########################################
function add_selected_cuts(model, θ, sys::system,tmax::Float64,P::Int64,constraints::Vector{Any},λ::Vector{Float64})
    success = attempts = 0 ; 
    max_violation = 0;
    while success<P && attempts<MAX_NUMBER_OF_ATTEMPTS
         t = Random.rand()*tmax;
         x = ((sys.xmax.-sys.xmin) .* [Random.rand() for j in 1:nx]) .+ sys.xmin;
         if is_feasible(sys,x)
            attempts+=1
            hmin,u = Hmin(sys,vcat(t,x),λ) ; 
            value = hmin + 1;
            if value < -1e-3
                max_violation = min(max_violation,value);
                append!(constraints, [add_hamiltonian_constraint(model,θ, sys,t,x,u)]);
                success+=1;
            end
        end
    end
    println("Number of attempts = ",attempts);
    println("Violation = ",max_violation);
    return max_violation
 end
 
 function add_selected_terminal_cuts(model, θ, sys::system,tmax::Float64,xT::Vector{Float64},P::Int64,constraints::Vector{Any},λ::Vector{Float64})
    success = attempts = 0 ; 
    max_violation = 0;
    while success<P && attempts<MAX_NUMBER_OF_ATTEMPTS_FINAL
        attempts+=1
         t = Random.rand()*tmax;
         u = random_vector_unit_ball(sys.nx);
         x = xT + RADIUS_BALL * u;
         value = v(vcat(t,x),λ);
         if value > 1e-3
            max_violation = max(max_violation,value);
            append!(constraints, [add_terminal_constraint(model, θ, t, x)]);
            success+=1;
         end
    end
    println("(Terminal) Number of attempts = ",attempts);
    println("(Terminal) Violation = ",max_violation);
    return max_violation
 end
    
    


##################################################### Main function ###########################################################
function dual_solving(sys,x0,xT,traj_heur,tmax,ϵ,μ)
    best_traj = traj_heur
    model = Model(Gurobi.Optimizer);
    set_optimizer_attribute(model, "Method", 1)
    @variable(model, θ[1:N]);
    constraints = [];
    add_hamiltonian_constraints_on_trajectory(model, θ, sys, traj_heur,constraints);
    add_hamiltonian_constraints_random(model, θ, sys,tmax,N_RANDOM_INIT,constraints);
    add_terminal_constraints(model, θ,xT,tmax,N_TERMINAL,constraints);
    #add_terminal_gradient_constraints(model, θ,sys,xT,tmax,N_TERMINAL,constraints);
    #add_terminal_sdp_constraints(model, θ,xT,tmax,200,constraints);

    if OBJECTIVE_WITH_INTEGRAL
        set_objective_with_integral(model,θ,x0,μ,tmax,0.1);
    else
        set_objective(model,θ,x0,μ);
    end
    λ = 0.0;
    max_violation = -Inf;
    while max_violation < -ϵ
        optimize!(model);
        λ = [value(θ[i]) for i in 1:N];
        success, traj,tmax_new = vmin_trajectory(sys,x0,xT,tmax,λ)
        if success && (tmax_new < tmax)
            tmax = tmax_new
            best_traj = traj
        end
        println("Time best feasible control = ",tmax);
        println("Value v(0,x0) = ",v(vcat(0.0,x0),λ));
        add_hamiltonian_constraints_on_trajectory(model, θ, sys, traj,constraints,λ);
        max_violation = add_selected_cuts(model, θ, sys,tmax,N_SELECTED_CONSTRAINTS,constraints,λ)

        if CIRCLE_FINAL_SET
            add_selected_terminal_cuts(model, θ, sys,tmax,xT,N_SELECTED_CONSTRAINTS_FINAL,constraints,λ)
        end
    end

    optimize!(model);
    λ = [value(θ[i]) for i in 1:N];
    success,traj,tmax_new = vmin_trajectory(sys,x0,xT,tmax,λ)
    if success && (tmax_new < tmax)
        tmax = tmax_new
        best_traj = traj
    end
    println("Time feasible control = ",tmax)
    println("Value v(0,x0) = ",v(vcat(0.0,x0),λ));
    return best_traj, tmax,λ
end