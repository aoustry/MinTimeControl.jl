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

function plot_hmin_traj(sys::system,λ::Vector{Float64},traj::Vector{Any})
    plot();
    plot!([el[1] for el in traj],  [Hmin(sys,el[1:1+sys.nx],λ)[1] for el in traj], xlabel = "Time", ylabel = "Hmin" );
    png("plots/"*sys.name*"hmintraj.png");
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

function circleShape(h,k,r)
    θ = LinRange(0,2*π,500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end
function arccircleShape(h,k,r,α)
    θ = LinRange(0,α,500)
    [h .+ r*sin.(θ);0.0], [k .+ r*cos.(θ);0.0]
end

function plot_water_flow(sys::system,traj_heur,best_traj)
        plot();
        amin,amax = sys.xmin[1], sys.xmax[1];
        bmin,bmax = sys.xmin[2], sys.xmax[2];
        a = amin: (amax-amin)/12 : amax;
        b = bmin: (bmax-bmin)/20 : bmax;
        resx, resy,resa,resb = [],[],[],[]
        for x in a
            for y in b
                append!(resa,x)
                append!(resb,y)
                append!(resx,0.17*sin(pi*y))
                append!(resy,0)
            end
        end
        quiver(resa,resb,quiver=(resx,resy),xlims = (amin,amax),ylims = (bmin,bmax),xlabel = "x₁",ylabel = "x₂",aspect_ratio=:equal);
        plot!(circleShape(0,-1,0.02), seriestype = [:shape,],lw = 0.5,c=:red,linecolor = :black, label = "Starting point", legend = :outertop,aspect_ratio=:equal)

        plot!(circleShape(0,0,TARGET_TOLERANCE), seriestype = [:shape,],lw = 0.5,c=:gray,linecolor = :black, label = "Target set K", legend = :outertop, fillalpha = 0.2,aspect_ratio=:equal)
        
        #plot!([el[2] for el in traj_heur], [el[3] for el in traj_heur],label = "Heuristic trajectory (t=1.261)", legend = :outertop,aspect_ratio=:equal );
        #plot!([el[2] for el in best_traj], [el[3] for el in best_traj],label = "Closed-loop trajectory (t=1.100)", legend = :outertop,xlabel = "x₁",ylabel = "x₂",aspect_ratio=:equal );
    
        xlims!(-0.8, 0.8)

        png("plots/waterflow.png");
end

function plot_control_zermelo(sys::system,traj_heur,best_traj)
    plot();
    
    plot!(circleShape(0,-1,0.02), seriestype = [:shape,],lw = 0.5,c=:red,linecolor = :black, label = "Starting point", legend = :outertop,aspect_ratio=:equal)

    plot!(circleShape(0,0,TARGET_TOLERANCE), seriestype = [:shape,],lw = 0.5,c=:gray,linecolor = :black, label = "Target set K", legend = :outertop, fillalpha = 0.2,aspect_ratio=:equal)
    plot!([el[2] for el in traj_heur], [el[3] for el in traj_heur],label = "Heuristic trajectory (t=1.261)", legend = :outertop,aspect_ratio=:equal );
    plot!([el[2] for el in best_traj], [el[3] for el in best_traj],label = "Closed-loop trajectory (t=1.100)", legend = :outertop,xlabel = "x₁",ylabel = "x₂",aspect_ratio=:equal );
    xlims!(-0.8, 0.8)
    png("plots/trajzermelo.png");
end

function plot_control_regatta(sys::system,traj_heur,best_traj)
    plot();
    plot!(circleShape(0,-1,0.02), seriestype = [:shape,],lw = 0.5,c=:red,linecolor = :black, label = "Starting point", legend = :outertop,aspect_ratio=:equal)
    plot!(circleShape(0,0,TARGET_TOLERANCE), seriestype = [:shape,],lw = 0.5,c=:gray,linecolor = :black, label = "Target set K", legend = :outertop, fillalpha = 0.2,aspect_ratio=:equal)
    plot!([el[2] for el in traj_heur], [el[3] for el in traj_heur],label = "Heuristic trajectory (t=1.278)", legend = :outertop,aspect_ratio=:equal );
    plot!([el[2] for el in best_traj], [el[3] for el in best_traj],label = "Closed-loop trajectory (t=0.913)", legend = :outertop,xlabel = "x₁",ylabel = "x₂",aspect_ratio=:equal );
    xlims!(-0.8, 0.8)
    png("plots/trajregatta.png");
end

function plotpolar(sys::system)
    plot()
    θ = LinRange(-π,π,500)

    r(α::AbstractFloat) = abs(sin(sys.polar_coef*α)) ;
    plot!([r(α) for α in θ] .*sin.(θ),  [r(α) for α in θ] .* cos.(θ),label = "Polar", legend = false,aspect_ratio=:equal );
    plot!([0.0,0.0],[0.0,.5],legend = false,color = :gray)
    θ = 2.0
    plot!([0.0,r(θ)*sin(θ)],[0.0,r(θ)*cos(θ)],legend = false,color = :gray)
    plot!(arccircleShape(0.0,0.0,0.15,θ), seriestype = [:shape,],lw = 0.5,c=:red,linecolor = :grey, legend = false,fillalpha = 0.2,aspect_ratio=:equal)

    png("plots/polar.png");
end

function plotwind(sys::system,t)
    plot();
    amin,amax = sys.xmin[1], sys.xmax[1];
    bmin,bmax = sys.xmin[2], sys.xmax[2];
    a = amin: (amax-amin)/8 : amax;
    b = bmin: (bmax-bmin)/8 : bmax;
    resx, resy,resa,resb = [],[],[],[]
    for x in a
        for y in b
            append!(resa,x)
            append!(resb,y)
            angle = 0.5*pi*(1-0.4*t)
            append!(resx,0.2*cos(angle+π))
            append!(resy,0.2*sin(angle+π))
        end
    end
    quiver(resa,resb,quiver=(resx,resy),xlims = (amin,amax),ylims = (bmin,bmax),xlabel = "x₁",ylabel = "x₂",aspect_ratio=:equal);
    #plot!(circleShape(0,-1,0.02), seriestype = [:shape,],lw = 0.5,c=:red,linecolor = :black, label = "Starting point", legend = :outertop,aspect_ratio=:equal)

    #plot!(circleShape(0,0,TARGET_TOLERANCE), seriestype = [:shape,],lw = 0.5,c=:gray,linecolor = :black, label = "Target set K", legend = :outertop, fillalpha = 0.2,aspect_ratio=:equal)
    
   #xlims!(-0.8, 0.8)

    png("plots/wind"*string(t)*".png");
end