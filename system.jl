using LinearAlgebra
using Optim
using Ipopt
using EAGO

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
        if (is_feasible(sys,x))==false
            return false, trajectory,t
        end 
        Δ = xT.-x;
        if norm(Δ) < TARGET_TOLERANCE
            return true,trajectory, t
        end
        u = heuristic_control(sys,t,x,Δ);
        append!(trajectory,[vcat(t,x,u)]);
        x = x .+ dt*f(sys,t,x,u);
        t += dt;
    end
    return false, trajectory, t
end

function vmin_trajectory(sys::system,x0::Vector{Float64},xT::Vector{Float64},Tmax::Float64,λ::Vector{Float64})
    t = 0.0;
    x = copy(x0);
    trajectory = [];
    while t <Tmax
        if  (is_feasible(sys,x))==false
            print(false,x)
            return false, trajectory,t
        end 
        if norm(xT.-x) < TARGET_TOLERANCE
            print(true,x)
            return true,trajectory, t
        elseif norm(xT.-x) < 1.2*TARGET_TOLERANCE
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
    print(false,x)
    return false,trajectory, t
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

function flow(sys::zermelo_boat,t::AbstractFloat,x::AbstractArray)
    @assert length(x)==sys.nx,
    return [-sys.flow_strength*sin(pi*x[2])*(1+sys.time_increasing_flow*t),0] ; 
end

function f(sys::zermelo_boat,t::AbstractFloat,x::AbstractArray,u::AbstractArray)
    return  u - flow(sys, t,x) ;
end

function argmin(sys::zermelo_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    return  -g/norm(g);
end

function heuristic_control(sys::zermelo_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    flow_vector = flow(sys,t,x);
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

function is_feasible(sys::zermelo_boat,x::Vector{Float64})
    return true
end

function certify_hjb(sys::zermelo_boat,λ::AbstractArray,tmax::Float64,optimizer)
    N = length(λ);
    model_cert = Model(optimizer);
    @variable(model_cert, 0<= t <= tmax);
    @variable(model_cert, sys.xmin[i] <= x[i = 1:2]<=sys.xmax[i]);
    @variable(model_cert, -π<= u <= π);
    #@NLconstraint(model_cert, u[1]^2 + u[2]^2==1.0);
    @NLexpression(model_cert, ft, 1.0);
    @NLexpression(model_cert, fx1, cos(u) + sys.flow_strength*sin(pi*x[2])*(1+sys.time_increasing_flow*t));
    @NLexpression(model_cert, fx2, sin(u));
    
    g1(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[1].α for i in 1:N)
    g2(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[2].α for i in 1:N)
    g3(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[3].α for i in 1:N)

    register(model_cert, :g1, 3, g1; autodiff = true)
    register(model_cert, :g2, 3, g2; autodiff = true)
    register(model_cert, :g3, 3, g3; autodiff = true)

    @NLobjective(model_cert, Min, 1.0 + g1(t,x[1],x[2])*ft + g2(t,x[1],x[2])*fx1 + g3(t,x[1],x[2])*fx2);
    optimize!(model_cert);
    return objective_value(model_cert),objective_bound(model_cert),[value(t),value(x[1]),value(x[2]),value(u)]
end

function certify_hjb_final(sys::zermelo_boat,λ::AbstractArray,tmax::Float64,xT::AbstractArray,optimizer)
    model_cert = Model(optimizer);
    @variable(model_cert, 0<= t <= tmax);
    @variable(model_cert, sys.xmin[i] <= x[i = 1:2]<=sys.xmax[i]);
    @NLconstraint(model_cert, (x[1]-xT[1])^2+(x[2]-xT[2])^2 <= TARGET_TOLERANCE^2);
    vstar(t::T, x1::T, x2:: T) where {T<:Real} = v([t,x1,x2],λ);
    register(model_cert, :vstar, 3, vstar; autodiff = true)
    @NLobjective(model_cert, Max, vstar(t,x[1],x[2]));
    optimize!(model_cert);
    return objective_value(model_cert),objective_bound(model_cert)
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

function wind_speed(sys::toy_boat,t::AbstractFloat,x::AbstractArray)
    return 2 + t
end

function wind_angle(sys::toy_boat,t::AbstractFloat,x::AbstractArray)
    return 0.5*pi*(1-0.4*t)
end

function polar(sys::toy_boat,rel_angle::AbstractFloat)
    return abs(sin(sys.polar_coef*rel_angle))
end


function f(sys::toy_boat,t::AbstractFloat,x::AbstractArray,u::AbstractArray)
    @assert abs(norm(u)-1)<1e-6;
    rel_angle = angle(u[1] + im*u[2]);
    wa = wind_angle(sys,t,x);
    #relative_angle = wind_angle(sys,t,x) - heading;
    #relative_angle = (relative_angle + pi)%(2*pi) - pi;
    r = wind_speed(sys,t,x)*polar(sys,rel_angle);
    return [r*cos(rel_angle+wa),r*sin(rel_angle+wa)];
end

#function fα1(sys::toy_boat,t::AbstractFloat,x::AbstractArray,α::AbstractFloat)
#    heading = α;
#    relative_angle = wind_angle(sys,t,x) - heading;
#    relative_angle = (relative_angle + pi)%(2*pi) - pi;
#    r = wind_speed(sys,t,x)*polar(sys,relative_angle);
#    return r*cos(α);
#end

function argmin(sys::toy_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    res = optimize(theta -> f(sys,t,x,[cos(theta),sin(theta)])'*g, -2*pi, 2*pi,Brent());
    θ = Optim.minimizer(res)[1];
    return [cos(θ),sin(θ)]
end


function heuristic_control(sys::toy_boat,t::Float64,x::Vector{Float64},g::Vector{Float64})
    target_angle = angle(g[1] + im*g[2]);
    wa = wind_angle(sys,t,x);
    θ = target_angle - wa;
    return [cos(θ),sin(θ)]
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

function is_feasible(sys::toy_boat,x::Vector{Float64})
    return true
end

function certify_hjb(sys::toy_boat,λ::AbstractArray,tmax::Float64,optimizer)
    N = length(λ);
    model_cert = Model(optimizer);
    @variable(model_cert, 0<= time <= tmax);
    @variable(model_cert, sys.xmin[i] <= x[i = 1:2]<=sys.xmax[i]);
    @variable(model_cert, -pi<= α <= pi);
    @variable(model_cert, -pi<= θ <= pi);
    ws(t::T) where {T<:Real} = 2 + t;
    polar(δ::T) where {T<:Real} =  (abs(sin(sys.polar_coef*δ)));
    register(model_cert, :ws, 1, ws; autodiff = true)
    register(model_cert, :polar, 1, polar; autodiff = true)
    g1(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[1].α for i in 1:N)
    g2(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[2].α for i in 1:N)
    g3(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[3].α for i in 1:N)
    register(model_cert, :g1, 3, g1; autodiff = true)
    register(model_cert, :g2, 3, g2; autodiff = true)
    register(model_cert, :g3, 3, g3; autodiff = true)
    @constraint(model_cert,α+0.5*pi*(1-0.4*time) == θ);
    @NLobjective(model_cert, Min, 1.0 + g1(time,x[1],x[2]) + g2(time,x[1],x[2])*ws(time)*cos(θ)*polar(α) 
                                                        + g3(time,x[1],x[2])*ws(time)*sin(θ)*polar(α))  ;
    optimize!(model_cert);
    return objective_value(model_cert),objective_bound(model_cert),[value(time),value(x[1]),value(x[2]),cos(value(α)),sin(value(α))]

end
 

#=  function certify_hjb(sys::toy_boat,λ::AbstractArray,tmax::Float64,optimizer)
    N = length(λ);
    model_cert = Model(optimizer);
    @variable(model_cert, 0<= time <= tmax);
    @variable(model_cert, sys.xmin[i] <= x[i = 1:2]<=sys.xmax[i]);
    @variable(model_cert, -pi<= α <= pi);
    ws(t::T) where {T<:Real} = 2 + t;
    wa(t::T) where {T<:Real} = 0.5*pi*(1-0.4*t);
    polar(δ::T) where {T<:Real} =  (abs(sin(sys.polar_coef*δ)));
    diff_angle(δ1::T,δ2::T) where {T<:Real} =  if (-π <= (δ1 - δ2)<=π) δ1 - δ2 elseif (-2*π <= (δ1 - δ2)<=-π) δ1 - δ2 + 2*π else  δ1 - δ2 - 2*π end  ;
    register(model_cert, :ws, 1, ws; autodiff = true)
    register(model_cert, :wa, 1, wa; autodiff = true)
    register(model_cert, :polar, 1, polar; autodiff = true)
    register(model_cert, :diff_angle, 2, diff_angle; autodiff = true)    
    g1(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[1].α for i in 1:N)
    g2(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[2].α for i in 1:N)
    g3(t::T, x1::T, x2:: T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2])[3].α for i in 1:N)
    register(model_cert, :g1, 3, g1; autodiff = true)
    register(model_cert, :g2, 3, g2; autodiff = true)
    register(model_cert, :g3, 3, g3; autodiff = true)
    @NLobjective(model_cert, Min, 1.0 + g1(time,x[1],x[2]) + g2(time,x[1],x[2])*ws(time)*cos(α)*polar(diff_angle(wa(time),α)) 
                                                        + g3(time,x[1],x[2])*ws(time)*sin(α)*polar(diff_angle(wa(time),α)))  ;
    optimize!(model_cert);
    return objective_value(model_cert),objective_bound(model_cert),[value(time),value(x[1]),value(x[2]),cos(value(α)),sin(value(α))]

end =#

function certify_hjb_final(sys::toy_boat,λ::AbstractArray,tmax::Float64,xT::AbstractArray,optimizer)
    model_cert = Model(optimizer);
    @variable(model_cert, 0<= t <= tmax);
    @variable(model_cert, sys.xmin[i] <= x[i = 1:2]<=sys.xmax[i]);
    @NLconstraint(model_cert, (x[1]-xT[1])^2+(x[2]-xT[2])^2 <= TARGET_TOLERANCE^2);
    vstar(t::T, x1::T, x2:: T) where {T<:Real} = v([t,x1,x2],λ);
    register(model_cert, :vstar, 3, vstar; autodiff = true)
    @NLobjective(model_cert, Max, vstar(t,x[1],x[2]));
    optimize!(model_cert);
    return objective_value(model_cert),objective_bound(model_cert)
end 

################################################ Brockett integrator test case  ###############################################################
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

function is_feasible(sys::brockett_integrator,x::Vector{Float64})
    return true
end

################################################ Generalized Brockett integrator test case  ###############################################################
struct gen_brockett_integrator <: system
    name::String
    nx::Integer;
    nu::Integer;
    xmax::Vector{Float64};
    umax::Vector{Float64};
    xmin::Vector{Float64};
    umin::Vector{Float64};
end



function f(sys::gen_brockett_integrator,t::Float64,x::Vector{Float64},u::Vector{Float64})
    vector = brockett_vector(sys,t,x);
    return  vcat(u,vector'*u) ;
end

function argmin(sys::gen_brockett_integrator,t::Float64,x::Vector{Float64},g::Vector{Float64})
    vector = brockett_vector(sys,t,x);
    v = g[1:sys.nu]+g[sys.nu+1]*vector;
    return  -v/norm(v);
end

function heuristic_control(sys::gen_brockett_integrator,t::Float64,x::Vector{Float64},g::Vector{Float64})
    if abs(xT[sys.nx]-x[sys.nx])>1e-2
        g = vcat(zeros(sys.nu),(xT[sys.nx]-x[sys.nx])/abs(xT[sys.nx]-x[sys.nx]));
    end
    vector = brockett_vector(sys,t,x);
    v = g[1:sys.nu]+g[sys.nu+1]*vector;
    return  v/norm(v);
end

function random_control(sys::gen_brockett_integrator,t::Float64,x::Vector{Float64})
    v = zeros(sys.nu)
    while norm(v) < .0001
        v = randn(sys.nu);
    end
    v = v / norm(v) 
end

function is_feasible(sys::gen_brockett_integrator,x::Vector{Float64})
    return true
end

function certify_hjb_final(sys::gen_brockett_integrator,λ::AbstractArray,tmax::Float64,xT::AbstractArray,optimizer)
    model_cert = Model(optimizer);
    @variable(model_cert, 0<= t <= tmax);
    @variable(model_cert, sys.xmin[i] <= x[i = 1:6]<=sys.xmax[i]);
    @NLconstraint(model_cert, (x[1]-xT[1])^2+(x[2]-xT[2])^2+(x[3]-xT[3])^2+(x[4]-xT[4])^2+(x[5]-xT[5])^2+(x[6]-xT[6])^2 <= TARGET_TOLERANCE^2);
    vstar(t::T, x1::T, x2:: T, x3::T,x4::T,x5::T,x6::T) where {T<:Real} = v([t,x1,x2,x3,x4,x5,x6],λ);
    register(model_cert, :vstar, 6, vstar; autodiff = true)
    @NLobjective(model_cert, Max, vstar(t,x[1],x[2],x[3],x[4],x[5],x[6]));
    optimize!(model_cert);
    return objective_value(model_cert),objective_bound(model_cert)
end


function certify_hjb(sys::gen_brockett_integrator,λ::AbstractArray,tmax::Float64,optimizer)
    N = length(λ);
    model_cert = Model(optimizer);
    @variable(model_cert, 0<= time <= tmax);
    @variable(model_cert, 0.99*sys.xmin[i]<= x[i = 1:sys.nx]<=0.99*sys.xmax[i]);
    @variable(model_cert, -1.0<= u[i=1:sys.nu] <= 1.0);
    @NLconstraint(model_cert, sum(u[i]*u[i] for i in 1:sys.nu)<=1.0);
    
    
    g1(t::T, x1::T, x2:: T,x3:: T,x4:: T,x5::T,x6::T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2,x3,x4,x5,x6])[1].α for i in 1:N)
    g2(t::T, x1::T, x2:: T,x3:: T,x4:: T,x5::T,x6::T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2,x3,x4,x5,x6])[2].α for i in 1:N)
    g3(t::T, x1::T, x2:: T,x3:: T,x4:: T,x5::T,x6::T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2,x3,x4,x5,x6])[3].α for i in 1:N)
    g4(t::T, x1::T, x2:: T,x3:: T,x4:: T,x5::T,x6::T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2,x3,x4,x5,x6])[4].α for i in 1:N)
    g5(t::T, x1::T, x2:: T,x3:: T,x4:: T,x5::T,x6::T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2,x3,x4,x5,x6])[5].α for i in 1:N)
    g6(t::T, x1::T, x2:: T,x3:: T,x4:: T,x5::T,x6::T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2,x3,x4,x5,x6])[6].α for i in 1:N)
    g7(t::T, x1::T, x2:: T,x3:: T,x4:: T,x5::T,x6::T) where {T<:Real} = sum(λ[i]*subs(symb∇[i],xvar=>[t,x1,x2,x3,x4,x5,x6])[7].α for i in 1:N)


    register(model_cert, :g1, 7, g1; autodiff = true)
    register(model_cert, :g2, 7, g2; autodiff = true)
    register(model_cert, :g3, 7, g3; autodiff = true)
    register(model_cert, :g4, 7, g4; autodiff = true)
    register(model_cert, :g5, 7, g5; autodiff = true)
    register(model_cert, :g6, 7, g6; autodiff = true)
    register(model_cert, :g7, 7, g7; autodiff = true)
    @NLobjective(model_cert, Min, 1.0 + g1(time,x[1],x[2],x[3],x[4],x[5],x[6]) + 
                                        g2(time,x[1],x[2],x[3],x[4],x[5],x[6]) * u[1] +
                                        g3(time,x[1],x[2],x[3],x[4],x[5],x[6]) * u[2] +
                                        g4(time,x[1],x[2],x[3],x[4],x[5],x[6]) * u[3] +
                                        g5(time,x[1],x[2],x[3],x[4],x[5],x[6]) * u[4] +
                                        g6(time,x[1],x[2],x[3],x[4],x[5],x[6]) * u[5] +  
                                         g7(time,x[1],x[2],x[3],x[4],x[5],x[6]) *(2*u[1]/(2+x[4])-x[1]*u[2]-cos(x[1]*x[3])*u[3]+exp(x[2])*u[4]+x[1]*x[2]*x[6]*u[5]));
   
    optimize!(model_cert);
    return objective_value(model_cert),objective_bound(model_cert),[value(time);[value(x[i]) for i in 1:sys.nx];[value(u[i]) for i in 1:sys.nu]]
end
