function [poly_rec, poly_rec_x] = fn_WENO3_rec_shw_trial_dambreak(xiq, x,dist_x_pl, w,flux,dflux)


%% flux splitting
a = max(abs(dflux(w)));  v=0.5*(flux(w)+a*w); u=circshift(0.5*(flux(w)-a*w),[0,-1]);
% We are seeking an approximation p(x_{iq}), of the function v(x_{iq}),
% where x_iq is in the interval [x_{i-1}, x_{i}].  We will look at two 
% candidate stencils: 
% {x_{i-2}, x_{i-1} x_{i}} and {x_{i-1}, x_{i}, x_{i+1}}
vm  = circshift(v,[0 1]);
vp  = circshift(v,[0 -1]);
vpp = circshift(v,[0,-2]);

xi_m  =  x - dist_x_pl;
xi    =  x;
xi_p  =  x + dist_x_pl;
xi_pp =  x + 2*dist_x_pl;

% Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}


p0n = vm .* (xiq - xi).* (xiq - xi_p)./(dist_x_pl.*(2*dist_x_pl)) - ...
      v .* (xiq - xi_m).*(xiq - xi_p)./((dist_x_pl).*(dist_x_pl)) +...
      vp .* (xiq - xi_m).*(xiq- xi)./((dist_x_pl).*(2* dist_x_pl));

p1n = v.* (xiq - xi_p).*(xiq - xi_pp)./((dist_x_pl).*(2*dist_x_pl))-...
      vp .* (xiq - xi).*(xiq - xi_pp)./((dist_x_pl).*(dist_x_pl))+...
      vpp .* (xiq - xi).*(xiq- xi_p)./((2*dist_x_pl).*(dist_x_pl));




% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];

             
v_arrs = [vm;v;vp;vpp];
v_ders = mat*v_arrs;

v_dm  =v_ders(1,:);
v_d  =v_ders(2,:);
v_dp  =v_ders(3,:);
v_dpp  =v_ders(4,:);
B0n = 4*(abs(v_dp - v_d) - abs(v_d-v_dm)).^2;
B1n = 4*(abs(v_dpp - v_dp) - abs(v_dp-v_d)).^2;

% linear weights
epsilon = 1E-6;
d0n = - (xiq - xi_pp)./(xi_pp - xi_m);

d1n = (xiq - xi_m)./(xi_pp - xi_m);

% Alpha weights 
pow = 1;
alpha0n = d0n./(epsilon + B0n).^pow;
alpha1n = d1n./(epsilon + B1n).^pow;
alphasumn = alpha0n + alpha1n;


% ENO stencils weights
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;

% Numerical value at 
poly_rec_n = w0n.*p0n + w1n.*p1n;


%% poly_rec_p

um  = circshift(u,[0 1]);
up  = circshift(u,[0 -1]);
upp = circshift(u,[0,-2]);

xi_m  =  x ;
xi    =  x + dist_x_pl;
xi_p  =  x + 2*dist_x_pl;
xi_pp =  x + 3*dist_x_pl;

% Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}


p0p = um .* (xiq - xi).* (xiq - xi_p)./(dist_x_pl.*(2*dist_x_pl)) - ...
      u .* (xiq - xi_m).*(xiq - xi_p)./((dist_x_pl).*(dist_x_pl)) +...
      up .* (xiq - xi_m).*(xiq- xi)./((dist_x_pl).*(2* dist_x_pl));

p1p = u.* (xiq - xi_p).*(xiq - xi_pp)./((dist_x_pl).*(2*dist_x_pl))-...
      up .* (xiq - xi).*(xiq - xi_pp)./((dist_x_pl).*(dist_x_pl))+...
      upp .* (xiq - xi).*(xiq- xi_p)./((2*dist_x_pl).*(dist_x_pl));




% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];

             
u_arrs = [um;u;up;upp];
u_ders = mat*u_arrs;

u_dm   = u_ders(1,:);
u_d    = u_ders(2,:);
u_dp   = u_ders(3,:);
u_dpp  = u_ders(4,:);
B0p = 4*(abs(u_dp - u_d) - abs(u_d-u_dm)).^2;
B1p = 4*(abs(u_dpp - u_dp) - abs(u_dp-u_d)).^2;

% linear weights
epsilon = 1E-6;
d0p = - (xiq - xi_pp)./(xi_pp - xi_m);

d1p = (xiq - xi_m)./(xi_pp - xi_m);

% Alpha weights 
pow = 1;
alpha0p = d0p./(epsilon + B0p).^pow;
alpha1p = d1p./(epsilon + B1p).^pow;
alphasump = alpha0p + alpha1p;

% ENO stencils weights
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;

% Numerical value at 
poly_rec_p = w0p.*p0p + w1p.*p1p;
% poly_rec = d0n.*p0n + d1n.*p1n;
% poly_rec_x = d0n.*p0n_x + d1n.*p1n_x + d0n_x.*p0n + d1n_x.*p1n;
poly_rec_x = (poly_rec_n-circshift(poly_rec_n,[0,1])+poly_rec_p-circshift(poly_rec_p,[0,1]))/dx;
poly_rec = poly_rec_n;
end
