using LinearAlgebra
using Optim

dt = 1e-4;
TARGET_TOLERANCE = 1e-2;

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
        end
        ∇ = ∇v(vcat(t,x),λ);
        u = argmin(sys,t,x,∇[2:end]);
        append!(trajectory,[vcat(t,x,u)]);
        x = x .+ dt*f(sys,t,x,u);
        t += dt;
    end
    return trajectory, t
end

################################################ Zermelo test case  ###############################################################
struct zermelo_boat <: system
    nx::Integer;
    nu::Integer;
    flow_strength::Float64;
    xmax::Vector{Float64};
    umax::Vector{Float64};
    xmin::Vector{Float64};
    umin::Vector{Float64};
end

function flow(sys::zermelo_boat,x::Vector{Float64})
    @assert length(x)==sys.nx,
    return [sys.flow_strength*sin(pi*(x[2]+0.5)),0] ; 
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
    nx::Integer
    nu::Integer
    polar_coef::Float64;
    xmax::Vector{Float64};
    umax::Vector{Float64};
    xmin::Vector{Float64};
    umin::Vector{Float64};
end

function wind_speed(sys::toy_boat,t::Float64,x::Vector{Float64})
    return 2.0
end

function wind_angle(sys::toy_boat,t::Float64,x::Vector{Float64})
    return 0.5*pi*(1-0.5*t)
end

function polar(sys::toy_boat,rel_angle::Float64)
    return abs(sin(sys.polar_coef*rel_angle))
end

#function aux(sys::toy_boat,g::Vector{Float64},winddir::Float64,theta::Float64)
#   diff_mod = (theta-winddir  + pi)%(2*pi) - pi
#    return (g[1]*cos(theta)+g[2]*sin(theta))*polar(sys,diff_mod)
#end

function f(sys::toy_boat,t::Float64,x::Vector{Float64},u::Vector{Float64})
    @assert abs(norm(u)-1)<1e-6;
    heading = angle(u[1] + im*u[2]);
    relative_angle = wind_angle(sys,t,x) - heading;
    relative_angle = (relative_angle + pi)%(2*pi) - pi;
    r = wind_speed(sys,t,x)*polar(sys,relative_angle);
    return r*u;
end

function argmin(sys::toy_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    windir = wind_angle(sys,t,x)
    x0 = 0;
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
