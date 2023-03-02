########################################### System definition #####################################################
function brockett_vector(sys::gen_brockett_integrator,t::Float64,x::Vector{Float64})
    #return [x[2],-x[1],sin(0.3*x[3]),exp(x[1]*x[2])]
    return [1/(1+x[4]),-x[1],-cos(x[1]*x[3]),exp(x[2])]
end
name_system = "gen_brockett"
nx, nu = 5,4;
xm,um,xM,uM = -1.0,-1.0,1.0,1.0
Tmax = 10.0;
x0 = [0.5,0.5,0.5,0.5,0.5]
xT = zeros(nx);
sys = gen_brockett_integrator(name_system,nx,nu,xM*ones(nx),uM*ones(nu),xm*ones(nx),um*ones(nu));
##################################################################################################################