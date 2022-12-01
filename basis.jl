using ForwardDiff
using DynamicPolynomials


################################################ Begin Basis parameterization ###############################################################
n = 1+nx;
degree = 8;
@polyvar xvar[1:n]
basis = monomials(xvar, 0:degree)
ϕ(y::AbstractArray) = [p(xvar=>y) for p in basis];
################################################ End Basis parameterization ###############################################################

N = length(basis)
∇ϕ(y::AbstractArray) = ForwardDiff.jacobian(ϕ,y);   
∇²ϕ(y::AbstractArray) = ForwardDiff.jacobian(x -> ∇ϕ(x),y);
g∇²ϕg(y::AbstractArray,g::Array) = ForwardDiff.derivative(t->∇ϕ(y+t*g)*g,0);
v(y::AbstractArray,λ::AbstractArray) = ϕ(y)'λ;
∇v(y::AbstractArray,λ::AbstractArray) = ForwardDiff.gradient(x -> v(x,λ),y);
∇²v(y::AbstractArray,λ::AbstractArray) = ForwardDiff.hessian(x -> v(x,λ),y);


