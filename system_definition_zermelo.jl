########################################### System definition #####################################################
param = 0.5
time_increasing_flow = 1.0
name_system = "zermelo_strength"*string(param)*"_increase"*string(time_increasing_flow);
nx, nu = 2,2;
xm,um,xM,uM = -1.0,-1.0,1.0,1.0
Tmax = 10.0;
x0 = [0.0,-1.0];
xT = [0,0.0]
sys = zermelo_boat(name_system,nx,nu,param,time_increasing_flow,[1.0,0.0],uM*ones(2),xm*ones(2),um*ones(2));
##################################################################################################################