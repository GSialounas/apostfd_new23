function [poly_rec, poly_rec_x] = fn_WENO3_rec_shw_nonu(xiq, x,dist_x_pl, w,flux,dflux)
% flux splitting
a = max(abs(dflux(w))); 
%pad w and xiq with ghost values
xiq = [xiq];%, xiq(end)+ (xiq(end)-xiq(end-1))];
w = [ w(2) w(1), w, w(end), w(end-1)];
dx_l = x(2)-x(1);
dx_r = x(end)-x(end-1);
x= [x(1) - 2*dx_l, x(1) - 1*dx_l, x, x(end) + dx_r, x(end) + 2*dx_r];
% the shift in the values of v and u (one to the left and one to the right)
% is to implement upwinding
v = w(1:end-1);
u = w(2:end);
a=max(abs(dflux(w))); v=0.5*(flux(v) + a*v);  u=0.5*(flux(u) - a*u); 
% v=0.5*(flux(w)+a*w); u=circshift(0.5*(flux(w)-a*w),[0,-1]);
% We are seeking an approximation p(x_{iq}), of the function v(x_{iq}),
% where x_iq is in the interval [x_{i-1}, x_{i}].  We will look at two 
% candidate stencils: 
% {x_{i-2}, x_{i-1} x_{i}} and {x_{i-1}, x_{i}, x_{i+1}}
vm = v(1:end-3);
vp = v(3:end-1);
vpp = v(4:end);
v= v(2:end-2);

x_downwind = x(1:end-1);
xi_m  =  x_downwind(1:end-3);
xi    =  x_downwind(2:end-2);
xi_p  =  x_downwind(3:end-1);
xi_pp = x_downwind(4:end);
% Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}
h_im = xi - xi_m;
h_i  = xi_p - xi;
h_ip = xi_pp - xi_p;

p0n = vm .* (xiq - xi).* (xiq - xi_p)./(h_im.*(h_im + h_i)) - ...
      v .* (xiq - xi_m).*(xiq - xi_p)./((h_im).*(h_i)) +...
      vp .* (xiq - xi_m).*(xiq- xi)./((h_im + h_i).*(h_i));

p1n = v.* (xiq - xi_p).*(xiq - xi_pp)./((h_i).*(h_i + h_ip))-...
      vp .* (xiq - xi).*(xiq - xi_pp)./((h_i).*(h_ip))+...
      vpp .* (xiq - xi).*(xiq- xi_p)./((h_i + h_ip).*(h_ip));


p0n_x =vm .* ((xiq - xi)+ (xiq - xi_p))./(h_im.*(h_im + h_i)) - ...
      v .* ((xiq - xi_m) + (xiq - xi_p))./((h_im).*(h_i)) +...
      vp .* ((xiq - xi_m)+ (xiq- xi))./((h_im + h_i).*(h_i));

p1n_x = v.* ((xiq - xi_p) + (xiq - xi_pp))./((h_i).*(h_i + h_ip))-...
      vp .* ((xiq - xi) + (xiq - xi_pp))./((h_i).*(h_ip))+...
      vpp .* ((xiq - xi) + (xiq- xi_p))./((h_i + h_ip).*(h_ip));

% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];
% 
%              
v_arrs = [vm;v;vp;vpp];
v_ders = mat*v_arrs;

v_dm  =v_ders(1,:);
v_d  =v_ders(2,:);
v_dp  =v_ders(3,:);
v_dpp  =v_ders(4,:);


B0n = 4*(abs(v_dp - v_d) - abs(v_d-v_dm)).^2;
B1n = 4*(abs(v_dpp - v_dp) - abs(v_dp-v_d)).^2;

% [nonu_ders, nonu_inds] =  fn_nonu_ders_non_periodic(x,vm, v, vp, vpp,h_im, h_i, h_ip);
% B0n = nonu_inds(1,:);
% B1n  =nonu_inds(2,:);
% linear weights
epsilon = 1E-6;
d0n = - (xiq - xi_pp)./(xi_pp - xi_m);
d1n = (xiq - xi_m)./(xi_pp - xi_m);

d0n_x = -1./(xi_pp - xi_m);
d1n_x = 1./(xi_pp - xi_m);

% Alpha weights 
pow = 1;
alpha0n = d0n./(epsilon + B0n).^pow;
alpha1n = d1n./(epsilon + B1n).^pow;
alphasumn = alpha0n + alpha1n;

alpha0n_x = d0n_x./(epsilon + B0n).^pow;
alpha1n_x = d1n_x./(epsilon + B1n).^pow;
alphasumn_x =  alpha0n_x +alpha1n_x;


% ENO stencils weights
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w0n_x = ((alpha0n_x).*(alphasumn) - alpha0n.*alphasumn_x )./ (alphasumn.^2);
w1n_x = ((alpha1n_x).*(alphasumn) - alpha1n.*alphasumn_x )./ (alphasumn.^2);


% Numerical value at 
poly_rec_n = w0n.*p0n + w1n.*p1n;
poly_rec_n_x = w0n.*p0n_x + w1n.*p1n_x + w0n_x.*p0n + w1n_x.*p1n;


%% poly_rec_p

um = u(1:end-3);
up = u(3:end-1);
upp = u(4:end);
u= u(2:end-2);

x_upwind = x(2:end);
xi_m  =  x_upwind(1:end-3);
xi    =  x_upwind(2:end-2);
xi_p  =  x_upwind(3:end-1);
xi_pp =  x_upwind(4:end);
% Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}
h_im = xi - xi_m;
h_i  = xi_p - xi;
h_ip = xi_pp - xi_p;

% um  = circshift(u,[0 1]);
% up  = circshift(u,[0 -1]);
% upp = circshift(u,[0,-2]);
% 
% xi_m  =  x ;
% xi    =  x + dist_x_pl;
% xi_p  =  x + 2*dist_x_pl;
% xi_pp =  x + 3*dist_x_pl;

% Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}


p0p = um .* (xiq - xi).* (xiq - xi_p)./(h_im.*(h_im + h_i)) - ...
      u .* (xiq - xi_m).*(xiq - xi_p)./((h_im).*(h_i)) +...
      up .* (xiq - xi_m).*(xiq- xi)./((h_im +h_i).*(h_i));

p1p = u.* (xiq - xi_p).*(xiq - xi_pp)./((h_i).*(h_i + h_ip))-...
      up .* (xiq - xi).*(xiq - xi_pp)./((h_i).*(h_ip))+...
      upp .* (xiq - xi).*(xiq- xi_p)./((h_i + h_ip).*(h_ip));

p0p_x = um .*( (xiq - xi) +  (xiq - xi_p))./(h_im.*(h_im + h_i)) - ...
      u .* ((xiq - xi_m) + (xiq - xi_p))./((h_im).*(h_i)) +...
      up .* ((xiq - xi_m) + (xiq- xi))./((h_im +h_i).*(h_i));

p1p_x = u.* ((xiq - xi_p) + (xiq - xi_pp))./((h_i).*(h_i + h_ip))-...
      up .* ((xiq - xi) + (xiq - xi_pp))./((h_i).*(h_ip))+...
      upp .* ((xiq - xi) + (xiq- xi_p)) ./((h_i + h_ip).*(h_ip));


% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];

%              
u_arrs = [um;u;up;upp];
u_ders = mat*u_arrs;

u_dm   = u_ders(1,:);
u_d    = u_ders(2,:);
u_dp   = u_ders(3,:);
u_dpp  = u_ders(4,:);
B0p = 4*(abs(u_dp - u_d) - abs(u_d-u_dm)).^2;
B1p = 4*(abs(u_dpp - u_dp) - abs(u_dp-u_d)).^2;
% % 
% [nonu_ders, nonu_inds] =  fn_nonu_ders_non_periodic(x,um, u, up, upp,h_im, h_i, h_ip);
% B0p = nonu_inds(1,:);
% B1p  =nonu_inds(2,:);
% % linear weights
epsilon = 1E-6;
d0p = - (xiq - xi_pp)./(xi_pp - xi_m);
d1p = (xiq - xi_m)./(xi_pp - xi_m);

d0p_x = -1./(xi_pp - xi_m);
d1p_x = 1./(xi_pp - xi_m);

% Alpha weights 
pow = 1;
alpha0p = d0p./(epsilon + B0p).^pow;
alpha1p = d1p./(epsilon + B1p).^pow;
alphasump = alpha0p + alpha1p;

alpha0p_x = d0p_x./(epsilon + B0p).^pow;
alpha1p_x = d1p_x./(epsilon + B1p).^pow;
alphasump_x =  alpha0p_x +alpha1p_x;

% ENO stencils weights
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;

w0p_x = ((alpha0p_x).*(alphasump) - alpha0p.*alphasump_x )./ (alphasump.^2);
w1p_x = ((alpha1p_x).*(alphasump) - alpha1p.*alphasump_x )./ (alphasump.^2);


% Numerical value at 
poly_rec_p = w0p.*p0p + w1p.*p1p;
poly_rec_p_x = w0p.*p0p_x + w1p.*p1p_x + w0p_x.*p0p + w1p_x.*p1p;% poly_rec_x = (poly_rec_n-circshift(poly_rec_n,[0,1])+poly_rec_p-circshift(poly_rec_p,[0,1]))/dx;

% poly_rec = d0n.*p0n + d1n.*p1n;
% poly_rec_x = d0n.*p0n_x + d1n.*p1n_x + d0n_x.*p0n + d1n_x.*p1n;
poly_rec_x =  (poly_rec_n_x+poly_rec_p_x);%(poly_rec_n(2:end)-poly_rec_n(1:end-1) + poly_rec_p(2:end)-poly_rec_p(1:end-1))/dx;
poly_rec = poly_rec_n;
end
