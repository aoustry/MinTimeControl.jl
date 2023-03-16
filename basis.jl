using ForwardDiff
using DynamicPolynomials
################################################ Begin Basis parameterization ###############################################################

#Polynomial basis
ϕ(y::AbstractArray) = [p(xvar=>y) for p in basis];
#Trigo basis
#ϕ(y::AbstractArray) = [ prod([cos(y[k]*(p.z[k]//2) + (p.z[k]%2) * pi/2) for k in 1:n]) for p in basis]
#N = length(basis)
#Combining?
#N = 2*length(basis)
#ϕ(y::AbstractArray) = vcat([p(xvar=>y) for p in basis], [ prod([cos(y[k]*(p.z[k]//2) + (p.z[k]%2) * pi/2) for k in 1:n]) for p in basis]);
################################################ End Basis parameterization ###############################################################
∇ϕ(y::AbstractArray) = ForwardDiff.jacobian(ϕ,y);  
symb∇ϕ(y::AbstractArray) = [subs(q,xvar=>y) for q in symb∇]
∇²ϕ(y::AbstractArray) = ForwardDiff.jacobian(x -> ∇ϕ(x),y);
g∇²ϕg(y::AbstractArray,g::Array) = ForwardDiff.derivative(t->∇ϕ(y+t*g)*g,0);
v(y::AbstractArray,λ::AbstractArray) = ϕ(y)'λ;
∇v(y::AbstractArray,λ::AbstractArray) = ForwardDiff.gradient(x -> v(x,λ),y);
symb∇v(y::AbstractArray,λ::AbstractArray) = λ'symb∇ϕ(y)
∇²v(y::AbstractArray,λ::AbstractArray) = ForwardDiff.hessian(x -> v(x,λ),y);


