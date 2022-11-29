using Plots

function plot_v(sys::system,λ::Vector{Float64},ub_t::Float64,name::String)
    for t in 0:ub_t/10:ub_t
        plot();
        amin,amax = sys.xmin[1], sys.xmax[1];
        bmin,bmax = sys.xmin[2], sys.xmax[2];
        a = amin: (amax-amin)/100 : amax;
        b = bmin: (bmax-bmin)/100 : bmax;
        contour!(a,b, (a,b)->v([t,a,b],λ));
        png("plots/"*name*string(t)*".png");
    end
end

function plot_v_terminal(sys::system,λ::Vector{Float64},x_T::Vector{Float64},tmax::Float64,name::String)
    plot();
    t = 0: (tmax)/100 : tmax;
    plot(t, (t)->v(vcat(t,xT),λ));
    png("plots/"*name*"_terminal.png");
end

function plot_∇v_terminal(sys::system,λ::Vector{Float64},x_T::Vector{Float64},tmax::Float64,name::String)
    plot();
    t = 0: (tmax)/100 : tmax;
    plot(t, (t)->norm(∇v(vcat(t,xT),λ)[2:sys.nx+1]));
    png("plots/"*name*"_grad_terminal.png");
end


function plot_traj(traj::Vector{Any},indices::Vector{Int},name::String)
    plot();
    plot!([el[indices[1]] for el in traj], [el[indices[2]] for el in traj] );
    png("plots/trajec"*name*".png");
end

function plot_v_traj(sys::system,λ::Vector{Float64},traj::Vector{Any},name::String)
    plot();
    plot!([el[1] for el in traj],  [v(el[1:1+sys.nx],λ) for el in traj], xlabel = "Time", ylabel = "Value function" );
    png("plots/"*name*".png");
end

