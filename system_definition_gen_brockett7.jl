########################################### System definition #####################################################
function brockett_vector(sys::gen_brockett_integrator,t::Float64,x::Vector{Float64})
    #return [x[2],-x[1],sin(0.3*x[3]),exp(x[1]*x[2]),x[2],x[1],x[3]^3,0,tan(x[7])]
    #return [x[1],-x[1],x[3],x[1],x[4],x[2]]
    return [x[3],x[4],x[5],x[6],x[1],x[2]]
end
name_system = "gen_brockett7"
nx, nu = 7,6;
xm,um,xM,uM = -1.0,-1.0,1.0,1.0
Tmax = 10.0;
x0 = [-0.2,0.4,0.2,0.2,0.2,-0.3,0.2]# [-0.3,0.2,-0.2] 
xT = zeros(nx);
sys = gen_brockett_integrator(name_system,nx,nu,xM*ones(nx),uM*ones(nu),xm*ones(nx),um*ones(nu));
##################################################################################################################