function [poly_rec, poly_rec_x] = fn_WENO3_rec(xiq, x,dist_x_pl, v)

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


% Spatial derivative of the polynomial 
p0n_x = vm .* ((xiq - xi) + (xiq - xi_p))./(dist_x_pl.*(2*dist_x_pl)) - ...
      v .* ((xiq - xi_m) + (xiq - xi_p))./((dist_x_pl).*(dist_x_pl)) +...
      vp .*( (xiq - xi_m) + (xiq- xi))./((dist_x_pl).*(2* dist_x_pl));
p1n_x = v.* ((xiq - xi_p) + (xiq - xi_pp))./((dist_x_pl).*(2*dist_x_pl))-...
      vp .* ((xiq - xi) + (xiq - xi_pp))./((dist_x_pl).*(dist_x_pl))+...
      vpp .*( (xiq - xi) + (xiq- xi_p))./((2*dist_x_pl).*(dist_x_pl));


% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
% smooth_mat = [-11, 18, -9, 2;...
%               -2 ,-3, 6, -1;...
%               1, -6,  3, 2 ;...
%               -2,9 ,-18, 11]
% B0n = (vm-v).^2;    
% B1n = (v-vp).^2;
% B0n = .25 * (abs(vp - vm) - abs(4* v - 3* vm - vp)).^2;
% B1n = .25 * (abs(vp - vm) - abs(4* v - vm - 3 * vp)).^2;
% 


dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];


             
v_arrs = [vm(1,:);v(1,:);vp(1,:);vpp(1,:)];
v_ders = mat*v_arrs;

v_dm = v_ders(1,:);
v_d  =v_ders(2,:);
v_dp  =v_ders(3,:);
v_dpp  =v_ders(4,:);

B0n = 4*(abs(v_dp - v_d) - abs(v_d-v_dm)).^2;
B1n = 4*(abs(v_dpp - v_dp) - abs(v_dp-v_d)).^2;
% 
% B0n = 4*(-vm +3*v -3*vp+vpp).^2;
% B1n =B0n;

% Constants
% d0n = 1/3; d1n = 2/3; epsilon = 1E-6;
epsilon = 1E-6;
d0n = - (xiq - xi_pp)./(xi_pp - xi_m);
d0n_x = -1./(xi_pp - xi_m);
d1n = (xiq - xi_m)./(xi_pp - xi_m);
d1n_x = 1./(xi_pp - xi_m);

% Alpha weights 
pow = 1;
alpha0n = d0n./(epsilon + B0n).^pow;
alpha1n = d1n./(epsilon + B1n).^pow;
alphasumn = alpha0n + alpha1n;

alpha0n_x = d0n_x./(epsilon + B0n).^pow;
alpha1n_x = d1n_x./(epsilon + B1n).^pow;
alphasumn_x = alpha0n_x + alpha1n_x;


% ENO stencils weights
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
% w0n_x = alpha0n_x./alphasumn;
% w1n_x = alpha1n_x./alphasumn;


w0n_x = ((alpha0n_x).*(alphasumn) - alpha0n.*alphasumn_x )./ (alphasumn.^2);
w1n_x = ((alpha1n_x).*(alphasumn) - alpha1n.*alphasumn_x )./ (alphasumn.^2);
% Numerical value at 
poly_rec = w0n.*p0n + w1n.*p1n;
poly_rec_x = w0n.*p0n_x + w1n.*p1n_x + w0n_x.*p0n + w1n_x.*p1n;


end
