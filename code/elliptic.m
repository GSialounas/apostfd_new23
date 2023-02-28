
%2nd order elliptic BVP -u_xx = f, u(-1) = u(1) = g, central difference
%approximation
clear all
close all

%Compute exact solution, derivative and RHS
syms x
uex = @(x) exp(-10*(x.^2));
Duex = matlabFunction( diff(uex(x)) );
f = matlabFunction( -diff(Duex(x)) );

for i = 1 : 15
    %Problem parameters
    N(i) = 2^(i+1);
    h = 2/(N(i)-1);
    nodes = linspace(-1,1,N(i));
    
    %System matrix
    A = sptoeplitz([2/h^2 -1/h^2 zeros(1, N(i)-2)]);
    
    %Enforce boundary conditions
    A(1,:) = 0.;
    A(1,1) = 1.;
    A(N(i),:) = 0.;
    A(N(i),N(i)) = 1.;
    b = f(nodes);
    b(1) = uex(nodes(1));
    b(N(i)) = uex(nodes(N(i)));
    
    %Solve system
    U = A\b';
    
    %Plot U as IU, the Lagrange interpolant.
    %plot(nodes,U)
    %hold on
    
    %Fit a cubic spline to the nodal values U
    xx = -1:.001:1;
    SU = spline(nodes,U);
    %plot(xx,ppval(SU,xx),'g')
    %plot(xx,uex(xx),'r')
    
    %Derivative of the spline
    [breaks,coefs,l,k,d] = unmkpp(SU);
    DSU = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
    
    %2nd derivative of the spline
    [breaks,coefs,l,k,d] = unmkpp(DSU);
    D2SU = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
    
    %residual for a posteriori estimate 
    err_func = @(x) (ppval(D2SU,x) + f(x)).^2;
    err(i) = sqrt(quad(err_func, -1,1));
    if i > 1
        EOC(i) = log(err_old/err(i))/log(2);
    end
    err_old = err(i);
end

loglog(N,err)
hold on
loglog(N,100*N.^(-2))