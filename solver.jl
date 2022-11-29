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
############################################################################################################

########################################## Individual constraints ##########################################
function add_hamiltonian_constraint(model,θ, sys::system,t::Float64,x::Vector{Float64},u::Vector{Float64})
    a =  ∇ϕ(vcat(t,x)) *vcat(1,f(sys,t,x,u));
    return @constraint(model,sum(a[i]*θ[i] for i in 1:N)>=-1);
end

function add_terminal_constraint(model, θ, t::Float64, xT::Vector{Float64})
    a = ϕ(vcat(t,xT));
    return @constraint(model,sum(a[i]*θ[i] for i in 1:N)<=0);
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

function add_terminal_sdp_constraints_random(model, θ,xT::Vector{Float64},tmax::Float64,P::Int64,constraints::Vector{Any})
    for t in tmax*0.5:tmax/(2*P):tmax
        for k in 1:10
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

function add_hamiltonian_constraints_on_trajectory(model, θ, sys::system, traj::Vector{Any},step::Int64,constraints::Vector{Any})
    for aux in 1:step:length(traj)
            vector = traj[aux]; 
            append!(constraints, [add_hamiltonian_constraint(model,θ, sys,vector[1],vector[2:nx+1],vector[nx+2:nx+nu+1])]);
    end
end

function add_hamiltonian_constraints_random(model, θ, sys::system,tmax::Float64,P::Int64,constraints::Vector{Any})
   for aux in 1:P
        t = Random.rand()*tmax;
        x = ((sys.xmax.-sys.xmin) .* [Random.rand() for j in 1:nx]) .+ sys.xmin;
        u = random_control(sys,t,x);
        append!(constraints, [add_hamiltonian_constraint(model,θ, sys,t,x,u)]);
   end
end
############################################################################################################




#ToBeDeleted:
#delete a constraint
#delete(model,constraint_by_name(model,"name"))