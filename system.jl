using LinearAlgebra

dt = 1e-4;
TARGET_TOLERANCE = 1e-2;

abstract type system end

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

#ToBeDeleted: ci-dessus: question de l'écart à -1 et de l'instantiation.

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
    xmax::Vector{Float64};
    umax::Vector{Float64};
    xmin::Vector{Float64};
    umin::Vector{Float64};
end