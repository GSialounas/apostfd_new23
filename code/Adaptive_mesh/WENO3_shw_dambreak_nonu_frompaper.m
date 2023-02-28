function [poly_rec, poly_rec_x] = WENO3_shw_dambreak_nonu_frompaper(w,flux,dflux,S,dx,eig_Jf_1,eig_Jf_2,F1,F2,x,t,L_domain)



xiq = x;
if length(w)==1
    a=max(abs(dflux(w)));
else
    eig_max_1 = max(abs(eig_Jf_1(w)));
    eig_max_2 = max(abs(eig_Jf_2(w)));
    a= max(eig_max_1,eig_max_2);
end
% xiq = x;
% if length(w)==1
%     a=max(abs(dflux(w)));
% else
%     eig_max_1 = max(abs(eig_Jf_1(w)));
%     eig_max_2 = max(abs(eig_Jf_2(w)));
%     a= max(eig_max_1,eig_max_2);
% end
% We are seeking an approximation p(x_{iq}), of the function v(x_{iq}),
% where x_iq is in the interval [x_{i-1}, x_{i}].  We will look at two 
% candidate stencils: 
% {x_{i-2}, x_{i-1} x_{i}} and {x_{i-1}, x_{i}, x_{i+1}}

dx_left = x(2) - x(1);
dx_right = x(end) - x(end-1);

w = [ w(:,2),w(:,1), w, w(:,end), w(:,end-1)];


x=  [ x(1) - 2*dx_left, x(1) - dx_left, x, x(end) + dx_right,x(end) + 2*dx_right];

% w = [ w(:,2),w(:,1), w, w(:,end)];
% x=  [ x(1)-2*dx_left,x(1) - dx_left, x, x(end) + dx_right];

% w = [ w(:,1), w, w(:,end), w(:,end-1),w(:,end-2)];
% x=  [  x(1) - dx_left, x, x(end) + dx_right,x(end) + 2*dx_right,x(end) + 3*dx_right];
% 
x_v= x(1:end-1);
x_u = x(2:end);

v = w(:,1:end-1);
u = w(:,2:end);
v=0.5*(flux(v) + a*v);  u=0.5*(flux(u) - a*u); 
% 
% x_u= x(1:end-1);
% x_v = x(2:end);
% 
% u = w(:,1:end-1);
% v = w(:,2:end);
% 
% x_u= x;
% x_v = x;
% 
% u = w;%(:,1:end-1);
% v = w;%(:,2:end);
% 
% u=0.5*(flux(u) + a*u);  v=0.5*(flux(v) - a*v); 
% vmm = circshift(v,[0,2]);
% vm  = circshift(v,[0 1]);
% vp  = circshift(v,[0 -1]);
% vpp = circshift(v,[0,-2]);


vm = v(:,1:end-3);
vp = v(:,3:end-1);
vpp = v(:,4:end);
v = v(:,2:end-2);

xi_m = x_v(1:end-3);
xi = x_v(2:end-2);
xi_p = x_v(3:end-1);
xi_pp = x_v(4:end);

h_im = xi - xi_m;
h_i = xi_p-xi;
h_ip = xi_pp-xi_p;

% xiq = xiq+h_i;


% h_im  =circshift(h_i,1);
% h_ip = circshift(h_i,-1);
% h_ipp = circshift(h_i,-2);

% h_imm = circshift(x,1) - circshift(x,2);
% h_imm(2)=L_domain-x(end);
% h_im  = circshift(x,0) - circshift(x,1);
% h_im(1) = L_domain-x(end);
% 
% h_i = circshift(x,-1) - circshift(x,0);
% h_i(end) =L_domain-x(end);%h_i(1);

% xi_m = x-h_im;%x - 2*dist_x_pl;
% xi    =  x;  
% xi_p  =  x+h_i;%x + dist_x_pl;
% xi_pp = x + h_i+ h_ip;

% Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}

p0n = vm.*((xiq - xi).*(xiq - xi_p))./(h_im.*(h_im+h_i)) - ...
      v.*((xiq-xi_m).*(xiq-xi_p))./(h_im.*h_i)+...
      vp.*((xiq - xi_m).*(xiq - xi))./((h_im + h_i).*(h_i));
  
p1n = v.*((xiq - xi_p).*(xiq - xi_pp))./(h_i.*(h_i+h_ip)) - ...
      vp.*((xiq-xi).*(xiq-xi_pp))./(h_i.*h_ip)+...
      vpp.*((xiq - xi).*(xiq - xi_p))./((h_i + h_ip).*(h_ip));

p0n_x = vm.*((xiq - xi)+(xiq - xi_p))./(h_im.*(h_im+h_i)) - ...
      v.*((xiq-xi_m)+(xiq-xi_p))./(h_im.*h_i)+...
      vp.*((xiq - xi_m)+(xiq - xi))./((h_im + h_i).*(h_i));
  
p1n_x = v.*((xiq - xi_p)+(xiq - xi_pp))./(h_i.*(h_i+h_ip)) - ...
      vp.*((xiq-xi)+(xiq-xi_pp))./(h_i.*h_ip)+...
      vpp.*((xiq - xi)+(xiq - xi_p))./((h_i + h_ip).*(h_ip));



% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];


             
v_arrs_h = [vm(1,:);v(1,:);vp(1,:);vpp(1,:)];
v_ders_h = mat*v_arrs_h;

v_arrs_hv = [vm(2,:);v(2,:);vp(2,:);vpp(2,:)];
v_ders_hv = mat*v_arrs_hv;

v_dm = [v_ders_h(1,:);v_ders_hv(1,:)];
v_d  =[v_ders_h(2,:);v_ders_hv(2,:)];
v_dp  =[v_ders_h(3,:);v_ders_hv(3,:)];
v_dpp  =[v_ders_h(4,:);v_ders_hv(4,:)];

[v_ders_h, v_inds_h] =  fn_nonu_ders_frompaper(h_im,h_i,h_ip,vm(1,:),v(1,:),vp(1,:),vpp(1,:));
[v_ders_hv, v_inds_hv] =  fn_nonu_ders_frompaper(h_im,h_i,h_ip,vm(2,:),v(2,:),vp(2,:),vpp(2,:));

B0n = [v_inds_h(1,:); v_inds_hv(1,:)];%4*(abs(v_dp - v_d) - abs(v_d-v_dm)).^2;
B1n = [v_inds_h(2,:); v_inds_hv(2,:)];%44*(abs(v_dpp - v_dp) - abs(v_dp-v_d)).^2;



% [nonu_ders, nonu_inds] =  fn_nonu_ders_nonu(x,v,L_domain);
% B0n  = nonu_inds(1,:);
% B1n = nonu_inds(2,:);
% % linear weights
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


%%%%
% upwind biased stencil


um = u(:,1:end-3);
up = u(:,3:end-1);
upp = u(:,4:end);
u = u(:,2:end-2);

xi_m = x_u(1:end-3);
xi = x_u(2:end-2);
xi_p = x_u(3:end-1);
xi_pp = x_u(4:end);

h_im = xi - xi_m;
h_i = xi_p-xi;
h_ip = xi_pp-xi_p;

% h_im  =circshift(h_i,1);
% h_ip = circshift(h_i,-1);
% h_ipp = circshift(h_i,-2);

% h_imm = circshift(x,1) - circshift(x,2);
% h_imm(2)=L_domain-x(end);
% h_im  = circshift(x,0) - circshift(x,1);
% h_im(1) = L_domain-x(end);
% 
% h_i = circshift(x,-1) - circshift(x,0);
% h_i(end) =L_domain-x(end);%h_i(1);

% xi_m = x-h_im;%x - 2*dist_x_pl;
% xi    =  x;  
% xi_p  =  x+h_i;%x + dist_x_pl;
% xi_pp = x + h_i+ h_ip;

% Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}
xiq = xiq ;

p0p = um.*((xiq - xi).*(xiq - xi_p))./(h_im.*(h_im+h_i)) - ...
      u.*((xiq-xi_m).*(xiq-xi_p))./(h_im.*h_i)+...
      up.*((xiq - xi_m).*(xiq - xi))./((h_im + h_i).*(h_i));
  
p1p = u.*((xiq - xi_p).*(xiq - xi_pp))./(h_i.*(h_i+h_ip)) - ...
      up.*((xiq-xi).*(xiq-xi_pp))./(h_i.*h_ip)+...
      upp.*((xiq - xi).*(xiq - xi_p))./((h_i + h_ip).*(h_ip));

p0p_x = um.*((xiq - xi)+(xiq - xi_p))./(h_im.*(h_im+h_i)) - ...
      u.*((xiq-xi_m)+(xiq-xi_p))./(h_im.*h_i)+...
      up.*((xiq - xi_m)+(xiq - xi))./((h_im + h_i).*(h_i));
  
p1p_x = u.*((xiq - xi_p)+(xiq - xi_pp))./(h_i.*(h_i+h_ip)) - ...
      up.*((xiq-xi)+(xiq-xi_pp))./(h_i.*h_ip)+...
      upp.*((xiq - xi)+(xiq - xi_p))./((h_i + h_ip).*(h_ip));



% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];


             
u_arrs_h = [um(1,:);u(1,:);up(1,:);upp(1,:)];
u_ders_h = mat*u_arrs_h;

u_arrs_hv = [um(2,:);u(2,:);up(2,:);upp(2,:)];
u_ders_hv = mat*u_arrs_hv;

u_dm   = [u_ders_h(1,:);u_ders_hv(1,:)];
u_d    = [u_ders_h(2,:);u_ders_hv(2,:)];
u_dp   = [u_ders_h(3,:);u_ders_hv(3,:)];
u_dpp  = [u_ders_h(4,:);u_ders_hv(4,:)];


[u_ders_h, u_inds_h] =  fn_nonu_ders_frompaper(h_im,h_i,h_ip,um(1,:),u(1,:),up(1,:),upp(1,:));
[u_ders_hv, u_inds_hv] =  fn_nonu_ders_frompaper(h_im,h_i,h_ip,um(2,:),u(2,:),up(2,:),upp(2,:));

B0p = [u_inds_h(1,:); u_inds_hv(1,:)];%4*(abs(v_dp - v_d) - abs(v_d-v_dm)).^2;
B1p = [u_inds_h(2,:); u_inds_hv(2,:)];%44*(abs(v_dpp - v_dp) - abs(v_dp-v_d)).^2;


% B0p = 4*(abs(u_dp - u_d) - abs(u_d - u_dm)).^2;
% B1p = 4*(abs(u_dpp - u_dp) - abs(u_dp - u_d)).^2;



% [nonu_ders, nonu_inds] =  fn_nonu_ders_nonu(x,v,L_domain);
% B0n  = nonu_inds(1,:);
% B1n = nonu_inds(2,:);
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
poly_rec_p_x = w0p.*p0p_x + w1p.*p1p_x + w0p_x.*p0p + w1p_x.*p1p;


poly_rec =  poly_rec_n + poly_rec_p;
poly_rec_x  =poly_rec_p_x +poly_rec_n_x;
% poly_rec = poly_rec(:,2:end);
% poly_rec_x = [poly_rec_x, poly_rec_x(:,end)];
% poly_rec_x  = poly_rec_n_x;
% poly_rec= [poly_rec, poly_rec(:,end) ];
% poly_rec_x = 1/dx_left*(poly_rec(:,2:end) - poly_rec(:,1:end-1));
end
