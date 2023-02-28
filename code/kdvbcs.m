% FD discretisation of u_t + uu_x + u_xxx = 0 on [xl,xr]
% BCs u(xl) = 0, u(xr) = 0, u_x(xr) = 0
% Note that the 3rd order term is treated implicitly with a midpoint
% timestepping scheme and the nonlinear term explicitly with an Euler
% timestepping scheme. This means the method is NOT unconditionally stable.
% There is a restriction on dt.

clear all
close all

showplot = 1; %boolean for graphics
T = 100.; % final time

% Some initial condition
ex = @(x,t) 0.5*(sech(0.5*((x-5)-t))).^2;

h = .05; % mesh size
% The spatial partition is [x0,x1,...xN,xN+1]. I take x0 = xl and xN =
% xr. xN+1 is a ghost point to help enforce the Neuman condition on the RH boundary.
x = [-20+h:h:20];

n = length(x);
dt = h^2; % timestep size

uold = ex(x,0)'; % set initial condition

if showplot
    fig = figure();
    lnh = plot(x,uold,'LineWidth',2); %approximate in blue
    
end
t = dt; %initialise time
i = 0;

% Dispersion matrix is the discretisation matrix associate to the 3rd order
% term
e = ones(n,1);
D = 1/(2*h^3)*spdiags([-e 2*e 0*e -2*e e], -2:2, n, n);

% To enforce the LH boundary use a one sided 2nd order difference quotient
% for the fist entry, this ensures u(xl) = 0
D(1,:) = 1/(2*h^3)*[10, -12, 6, -1, zeros(1,n-4)];

% For the RH boundary dirichlet term, add an entry
D(end,:) = 1/(2*h^3)*[zeros(1,n-3), -1, 2, 1];

% The nonlinear term is going to be treated explicitly. We only need an
% advection matrix. Use a 2nd order central approximation 
N = 1/(2*h)*spdiags([-1*e 0*e 1*e], -1:1, n, n);

while t < T+eps
    i = i + 1;

    u = (eye(n) + dt/2*D) \ (uold -dt/2*N*(uold.^2) - dt/2*D*uold);
        
    if showplot && mod(i,1) == 0
        set(lnh,'XData',x,'YData',u)
        drawnow
    end
    
    t = t + dt; %move in time
    uold = u;
    
end
