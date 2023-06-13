using JuMP
using SumOfSquares
using DynamicPolynomials
using CSDP
using LinearAlgebra # Needed for `dot`
include("basis.jl")
include("system.jl")
system_file_name = "system_definition_zermelo.jl";
include(system_file_name)

RADIUS_BALL = TARGET_TOLERANCE = 5*1e-2; #1e-2
dt = 1e-3; #1e-3




function vmin_trajectory_sos(sys::system,x0::Vector{Float64},xT::Vector{Float64},Tmax::Float64)
    t_elapsed = 0.0;
    x = copy(x0);
    trajectory = [];
    while t_elapsed <Tmax
        if  (is_feasible(sys,x))==false
            return false, trajectory,t_elapsed
        end 
        if norm(xT.-x) < TARGET_TOLERANCE
            return true,trajectory, t_elapsed
        elseif norm(xT.-x) < 1.2*TARGET_TOLERANCE
            Δ = xT.-x;
            u = heuristic_control(sys,t_elapsed,x,Δ);
        else
            ∇ = ∇valV(vcat(t_elapsed,x));
            u = argmin(sys,t_elapsed,x,∇[2:end]);
        end
        append!(trajectory,[vcat(t_elapsed,x,u)]);
        x = x .+ dt*f(sys,t_elapsed,x,u);
        t_elapsed += dt;
    end
    return false,trajectory, t_elapsed
end

for degs in [2,4,6,8]
time_sos = @elapsed begin    

nx = 4;
@polyvar t
@polyvar x[1:nx]; #(x3 = cos (pi*x2), x4 = sin(pi*x2))
@polyvar u[1:2]
valV(y::AbstractArray) = value(v_poly(t=>y[1],x[1]=>y[2],x[2]=>y[3],x[3]=>cos(pi*y[3]),x[4]=>sin(pi*y[3])));
∇valV(y::AbstractArray) = ForwardDiff.gradient(z -> valV(z),y);
monos = monomials([t;x], degs)
M = length(monos) + 1;
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);
model = SOSModel(solver);
@variable(model,λvar[1:M]);
v_poly = λvar[1]+polynomial(λvar[2:end], monos);

t0 = 0;
F = [u[1] + 0.5*(1+t)*x[4] ; u[2]; -pi*x[4]*u[2]; pi*x[3]*u[2]]; 
dVdt = differentiate(v_poly, t) + dot(differentiate(v_poly, x), F);

#Lyapunov constraint
ItimesXtimesU = @set x[1]>=-1 && x[1] <= 1 && x[2] <= 0 && x[2]>=-1 && x[3]>=-1 && x[3] <= 1 && x[4] <= 0 && x[4]>=-1 && u[1]^2 + u[2]^2 <=1 && t >= 0 && t<=1.3 && x[3]^2 + x[4]^2 == 1 ;
@constraint(model, Lyapunov_constraint, dVdt + 1 >= 0, domain = ItimesXtimesU)


#Terminal constraint
ItimesK = @set x[1]^2 + x[2]^2 <= 0.05^2 && x[3]>=-1 && x[3] <= 1 && x[4] <= 0 && x[4]>=-1 && t >= 0 && t<=1.3 && x[3]^2 + x[4]^2 == 1 ;
@constraint(model, terminal_constraint, v_poly <= 0, domain = ItimesK)

@objective(model, Max, v_poly(t=>0,x[1]=>0,x[2]=>-1,x[3]=>-1,x[4]=>0))
optimize!(model)

lb= objective_value(model);
λvalue = value.(λvar);
bool, traj, ub = vmin_trajectory_sos(sys,[0.0,-1.0],[0.0,0.0],1.3);

end

print(degs, " & ",lb," & ",ub ," & ",time_sos,"\\")

end