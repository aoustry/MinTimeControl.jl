using LinearAlgebra
using Optim



abstract type system end


function H(sys::system,t::Float64,x::Vector{Float64},u::Vector{Float64},λ::Vector{Float64})
    ∇ = ∇v(vcat(t,x),λ);
    return  ∇' *vcat(1,f(sys,t,x,u));
end

function Hmin(sys::system,y::Vector{Float64},λ::Vector{Float64})
    ∇ = ∇v(y,λ);
    u = argmin(sys,y[1],y[2:end],∇[2:end]);
    return  ∇' *vcat(1,f(sys,y[1],y[2:end],u)),u;
end


function heuristic_trajectory(sys::system,x0::Vector{Float64},xT::Vector{Float64},Tmax::Float64)
    t = 0.0;
    x = copy(x0);
    trajectory = [];
    while t <Tmax
        Δ = xT.-x;
        if norm(Δ) < TARGET_TOLERANCE
            return trajectory, t
        end
        u = heuristic_control(sys,t,x,Δ);
        append!(trajectory,[vcat(t,x,u)]);
        x = x .+ dt*f(sys,t,x,u);
        t += dt;
    end
    return trajectory, t
end

function vmin_trajectory(sys::system,x0::Vector{Float64},xT::Vector{Float64},Tmax::Float64,λ::Vector{Float64})
    t = 0.0;
    x = copy(x0);
    trajectory = [];
    while t <Tmax
        if norm(xT.-x) < TARGET_TOLERANCE
            return trajectory, t
        elseif norm(xT.-x) < 3*TARGET_TOLERANCE
            Δ = xT.-x;
            u = heuristic_control(sys,t,x,Δ);
        else
            ∇ = ∇v(vcat(t,x),λ);
            u = argmin(sys,t,x,∇[2:end]);
        end
        append!(trajectory,[vcat(t,x,u)]);
        x = x .+ dt*f(sys,t,x,u);
        t += dt;
    end
    return trajectory, t
end

################################################ Zermelo test case  ###############################################################
struct zermelo_boat <: system
    name::String
    nx::Integer;
    nu::Integer;
    flow_strength::Float64;
    time_increasing_flow::Float64;
    xmax::Vector{Float64};
    umax::Vector{Float64};
    xmin::Vector{Float64};
    umin::Vector{Float64};
end

function flow(sys::zermelo_boat,x::Vector{Float64})
    @assert length(x)==sys.nx,
    return [sys.flow_strength*sin(pi*(x[2]+0.5))*(1+sys.time_increasing_flow*x[1]),0] ; 
end

function f(sys::zermelo_boat,t::Float64,x::Vector{Float64},u::Vector{Float64})
    return  u - flow(sys, x) ;
end

function argmin(sys::zermelo_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    return  -g/norm(g);
end

function heuristic_control(sys::zermelo_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    flow_vector = flow(sys,x);
    a = norm(g)^2;
    b = 2*flow_vector'*g;
    c = norm(flow_vector)^2 - 1;
    disc = b^2-4*a*c;
    l = (-b+sqrt(disc))/(2*a);
    u = l*g+flow_vector;
    @assert abs(norm(u)-1)<1e-6;
    return u
end

function random_control(sys::zermelo_boat,t::Float64,x::Vector{Float64})
    angle = Random.rand() * 2* pi;
    return [cos(angle),sin(angle)]
end


################################################ toy_boat test case  ###############################################################
struct toy_boat <: system
    name::String
    nx::Integer
    nu::Integer
    polar_coef::Float64;
    xmax::Vector{Float64};
    umax::Vector{Float64};
    xmin::Vector{Float64};
    umin::Vector{Float64};
end

function wind_speed(sys::toy_boat,t::Float64,x::Vector{Float64})
    return 2.5
end

function wind_angle(sys::toy_boat,t::Float64,x::Vector{Float64})
    return 0.5*pi*(1-0.4*t)
end

function polar(sys::toy_boat,rel_angle::Float64)
    return abs(sin(sys.polar_coef*rel_angle))
end


function f(sys::toy_boat,t::Float64,x::Vector{Float64},u::Vector{Float64})
    @assert abs(norm(u)-1)<1e-6;
    heading = angle(u[1] + im*u[2]);
    relative_angle = wind_angle(sys,t,x) - heading;
    relative_angle = (relative_angle + pi)%(2*pi) - pi;
    r = wind_speed(sys,t,x)*polar(sys,relative_angle);
    return r*u;
end

function argmin(sys::toy_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    res = optimize(theta -> f(sys,t,x,[cos(theta),sin(theta)])'*g, -2*pi, 2*pi,Brent());
    θ = Optim.minimizer(res)[1];
    return [cos(θ),sin(θ)]
end
    

function heuristic_control(sys::toy_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    return g/norm(g)
end

function random_control(sys::toy_boat,t::Float64,x::Vector{Float64})
    angle = Random.rand() * 2* pi;
    return [cos(angle),sin(angle)]
end

function plot_test(sys::system,t::Float64)
    θ = -pi:0.01:pi;
    Y1 = [f(sys,t,[0.0,0.0],[cos(a),sin(a)])[1] for a in θ]
    Y2 = [f(sys,t,[0.0,0.0],[cos(a),sin(a)])[2] for a in θ]
    plot(Y1,Y2,aspect_ratio = 1);
    png("tests/"*sys.name*string(t)*".png")
end

function plot_test_direction(sys::system,t::Float64,g ::Vector{Float64})
    θ = -pi:0.01:pi;
    Y1 = [f(sys,t,[0.0,0.0],[cos(a),sin(a)])[1] for a in θ]
    Y2 = [f(sys,t,[0.0,0.0],[cos(a),sin(a)])[2] for a in θ]
    u = argmin(sys,t,[0.0,0.0],g)
    plot(Y1,Y2,aspect_ratio = 1);
    plot!([0.0,u[1]],[0.0,u[2]],aspect_ratio = 1)
    png("tests/"*sys.name*string(t)*".png")
end


################################################ Zermelo test case  ###############################################################
struct brockett_integrator <: system
    name::String
    nx::Integer;
    nu::Integer;
    xmax::Vector{Float64};
    umax::Vector{Float64};
    xmin::Vector{Float64};
    umin::Vector{Float64};
end



function f(sys::brockett_integrator,t::Float64,x::Vector{Float64},u::Vector{Float64})
    return  [u[1],u[2],u[1]*x[2]-u[2]*x[1]] ;
end

function argmin(sys::brockett_integrator,t::Float64,x::Vector{Float64},g::Vector{Float64})
    v = [g[1]+g[3]*x[2], g[2] - g[3]*x[1]];
    return  -v/norm(v);
end

function heuristic_control(sys::brockett_integrator,t::Float64,x::Vector{Float64},g::Vector{Float64})
    if abs(xT[3]-x[3])>1e-3
        g = [0.0,0.0,(xT[3]-x[3])/abs(xT[3]-x[3])]
    end
    v = [g[1]+g[3]*x[2], g[2] - g[3]*x[1]];
    return  v/norm(v);
end

function random_control(sys::brockett_integrator,t::Float64,x::Vector{Float64})
    angle = Random.rand() * 2* pi;
    return [cos(angle),sin(angle)]
end
