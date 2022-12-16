using Plots

function plot_v(sys::system,λ::Vector{Float64},ub_t::Float64)
    for t in 0:ub_t/10:ub_t
        plot();
        amin,amax = sys.xmin[1], sys.xmax[1];
        bmin,bmax = sys.xmin[2], sys.xmax[2];
        a = amin: (amax-amin)/100 : amax;
        b = bmin: (bmax-bmin)/100 : bmax;
        contour!(a,b, (a,b)->v([t,a,b],λ));
        png("plots/"*sys.name*string(t)*".png");
    end
end

function plot_v_terminal(sys::system,λ::Vector{Float64},x_T::Vector{Float64},tmax::Float64)
    plot();
    t = 0: (tmax)/100 : tmax;
    plot(t, (t)->v(vcat(t,xT),λ));
    png("plots/"*sys.name*"_terminal.png");
end

function plot_∇v_terminal(sys::system,λ::Vector{Float64},x_T::Vector{Float64},tmax::Float64)
    plot();
    t = 0: (tmax)/100 : tmax;
    plot(t, (t)->norm(∇v(vcat(t,xT),λ)[2:sys.nx+1]));
    png("plots/"*sys.name*"_grad_terminal.png");
end


function plot_traj(traj::Vector{Any},indices::Vector{Int})
    plot();
    plot!([el[indices[1]] for el in traj], [el[indices[2]] for el in traj] );
    png("plots/trajec"*sys.name*".png");
end

function plot_v_traj(sys::system,λ::Vector{Float64},traj::Vector{Any})
    plot();
    plot!([el[1] for el in traj],  [v(el[1:1+sys.nx],λ) for el in traj], xlabel = "Time", ylabel = "Value function" );
    png("plots/"*sys.name*".png");
end

function plot_hmin(sys::system,λ::Vector{Float64},ub_t::Float64)
    for t in 0:ub_t/20:ub_t
        plot();
        amin,amax = sys.xmin[1], sys.xmax[1];
        bmin,bmax = sys.xmin[2], sys.xmax[2];
        a = amin: (amax-amin)/100 : amax;
        b = bmin: (bmax-bmin)/100 : bmax;
        resx, resy = [],[]
        for x in a
            for y in b
                h,u = Hmin(sys,[t,x,y],λ)
                if h < -0.999
                    append!(resx,x)
                    append!(resy,y)
                end
            end
        end
        scatter(resx,resy,xlims = (amin,amax),ylims = (bmin,bmax));
        png("plots/"*sys.name*string(t)*"hmin.png");
    end
end

function plot_field(sys::system,λ::Vector{Float64},ub_t::Float64)
    counter = 0
    for t in 0:ub_t/20:ub_t
        plot();
        amin,amax = sys.xmin[1], sys.xmax[1];
        bmin,bmax = sys.xmin[2], sys.xmax[2];
        a = amin: (amax-amin)/30 : amax;
        b = bmin: (bmax-bmin)/30 : bmax;
        resx, resy,resa,resb = [],[],[],[]
        for x in a
            for y in b
                h,u = Hmin(sys,[t,x,y],λ)
                append!(resa,x)
                append!(resb,y)
                append!(resx,0.1*u[1])
                append!(resy,0.1*u[2])
            end
        end
        quiver(resa,resb,quiver=(resx,resy),xlims = (amin,amax),ylims = (bmin,bmax));
        png("plots/field"*sys.name*"_"*string(counter)*"_"*string(t)*"hmin.png");
        counter+=1;
    end
end