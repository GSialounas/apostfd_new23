function [poly_rec, poly_rec_x] = WENO3_shw_dambreak_nonu(w,flux,dflux,S,dx,eig_Jf_1,eig_Jf_2,F1,F2,x,t,L_domain)


%% flux splitting
% a = max(abs(dflux(w)));  v=0.5*(flux(w)+a*w); u=circshift(0.5*(flux(w)-a*w),[0,-1]);
% We are seeking an approximation p(x_{iq}), of the function v(x_{iq}),
% where x_iq is in the interval [x_{i-1}, x_{i}].  We will look at two 
% candidate stencils: 
% {x_{i-2}, x_{i-1} x_{i}} and {x_{i-1}, x_{i}, x_{i+1}}

xiq = x;
if length(w)==1
    a=max(abs(dflux(w)));
else
    eig_max_1 = max(abs(eig_Jf_1(w)));
    eig_max_2 = max(abs(eig_Jf_2(w)));
    a= max(eig_max_1,eig_max_2);
end
dx_left = x(2) - x(1);
dx_right = x(end) - x(end-1);
% 
w = [ w(:,2) w(:,1), w, w(:,end), w(:,end-1)];
x=  [x(1) - 2*dx_left, x(1) - dx_left, x, x(end) + dx_right,x(end) + 2*dx_right];
% w = [  w(:,1), w, w(:,end), w(:,end-1), w(:,end-2)];
% x=  [ x(1) - dx_left, x, x(end) + dx_right,x(end) + 2*dx_right,x(end) + 3*dx_right];

% w = [w(:,3), w(:,2) w(:,1), w, w(:,end)];
% x=  [x(1) - 3*dx_left, x(1) - 2*dx_left, x(1) - dx_left, x, x(end) + dx_right];

v = w(:,1:end-1);
v_ders = v; % this will be used for the smoothness indicator
x_v = x(1:end-1);

u = w(:,2:end);
u_ders = u; % this will be used for the smoothness indicator

x_u = x(2:end);

% v=0.5*(flux(w)+a*w); u=circshift(0.5*(flux(w)-a*w),[0,-1]);
v=0.5*(flux(v) + a*v);  u=0.5*(flux(u) - a*u); 


vmm = v(:,1:end-3);
vm = v(:,2:end-2);
vp = v(:,4:end);
v = v(:,3:end-1);

xi_mm = x_v(1:end-3);
xi_m = x_v(2:end-2);
xi_p = x_v(4:end);
xi = x_v(3:end-1);
xi_v = xi;

% 
% 
% vmm = circshift(v,[0,2]);
% vm  = circshift(v,[0 1]);
% vp  = circshift(v,[0 -1]);
% vpp = circshift(v,[0,-2]);

h_i = xi_p -xi;
h_im = xi - xi_m;
h_imm = xi_m -xi_mm;

% h_i = circshift(x,-1) - circshift(x,0);
% h_i(end) =L_domain-x(end);%h_i(1);
% h_im  = circshift(h_i,1);
% h_imm = circshift(h_i,2);
% h_ip = circshift(h_i,-1);

% h_imm = circshift(x,1) - circshift(x,2);
% h_imm(2)=L_domain-x(end);
% h_im  = circshift(x,0) - circshift(x,1);
% h_im(1) = L_domain-x(end);
% 
% h_i = circshift(x,-1) - circshift(x,0);
% h_i(end) =L_domain-x(end);%h_i(1);

% xi_mm = x-h_im-h_imm;%x - 2*dist_x_pl;
% xi_m  =  x-h_im;%x - dist_x_pl;
% xi    =  x;  
% xi_p  =  x+h_i;%x + dist_x_pl;

% Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}

p0n = vmm.*(xiq -xi_m).*(xiq- xi)./((h_imm).*(h_imm+h_im) )- ...
      vm .* (xiq -xi_mm).*(xiq- xi)./((h_imm).*(h_im) )+ ...
      v .* (xiq -xi_mm).*(xiq- xi_m)./((h_imm+h_im).*(h_im) );

p1n = vm.*(xiq - xi).*(xiq -xi_p)./((h_im).*(h_im+h_i))- ...
      v .*(xiq - xi_m).*(xiq - xi_p)./((h_im).*(h_i))+ ...
      vp.*(xiq - xi_m).*(xiq - xi)./((h_im+h_i).*(h_i));
  
% p0n = vm .* (xiq - xi).* (xiq - xi_p)./(dist_x_pl.*(2*dist_x_pl)) - ...
%       v .* (xiq - xi_m).*(xiq - xi_p)./((dist_x_pl).*(dist_x_pl)) +...
%       vp .* (xiq - xi_m).*(xiq- xi)./((dist_x_pl).*(2* dist_x_pl));
% 
% p1n = v.* (xiq - xi_p).*(xiq - xi_pp)./((dist_x_pl).*(2*dist_x_pl))-...
%       vp .* (xiq - xi).*(xiq - xi_pp)./((dist_x_pl).*(dist_x_pl))+...
%       vpp .* (xiq - xi).*(xiq- xi_p)./((2*dist_x_pl).*(dist_x_pl));
% 
p0n_x = vmm.*((xiq -xi_m) + (xiq- xi))./((h_imm).*(h_imm+h_im) )- ...
        vm .* ((xiq -xi_mm) + (xiq- xi))./((h_imm).*(h_im) )+ ...
        v .* ((xiq -xi_mm) + (xiq- xi_m))./((h_imm+h_im).*(h_im) );

p1n_x = vm.*((xiq - xi) + (xiq -xi_p))./((h_im).*(h_im+h_i))-...
      v .*((xiq - xi_m) + (xiq - xi_p))./((h_im).*(h_i))+ ...
      vp.*((xiq - xi_m) + (xiq - xi))./((h_im+h_i).*(h_i));
  


% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];

             
v_arrs_h = [vmm(1,:);vm(1,:);v(1,:);vp(1,:)];
v_ders_h = mat*v_arrs_h;

v_arrs_hv = [vmm(2,:);vm(2,:);v(2,:);vp(2,:)];
v_ders_hv = mat*v_arrs_hv;

v_dmm = [v_ders_h(1,:);v_ders_hv(1,:)];
v_dm = [v_ders_h(2,:);v_ders_hv(2,:)];
v_d  =[v_ders_h(3,:);v_ders_hv(3,:)];
v_dp  =[v_ders_h(4,:);v_ders_hv(4,:)];

% B0n = (vm-v).^2; 
% B1n = (v-vp).^2;

% B0n  = .25 * (abs(v - vmm) - abs(4*vm - 3*vmm +v)).^2;
% B1n = .25*(abs(v - vmm) - abs(4*vm - vmm -3*v)).^2;

B0n = 4*(abs(v_d - v_dm) - abs(v_dm-v_dmm)).^2;
B1n = 4*(abs(v_dp - v_d) - abs(v_d-v_dm)).^2;
% 
% [nonu_ders, nonu_inds] =  fn_nonu_ders_nonu_non_periodic(x_v,v_ders,L_domain);
% B0n  = nonu_inds(1,:);
% B1n = nonu_inds(2,:);
% linear weights
epsilon =1e-6;
d0n = - (xiq - xi_p)./(xi_p - xi_mm);

d1n = (xiq - xi_mm)./(xi_p - xi_mm);
d0n_x = -1./(xi_p - xi_mm);
d1n_x = 1./(xi_p - xi_mm);

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
% umm um u up and upp correspond to 
% vm v 

umm = u(:,1:end-3);
um = u(:,2:end-2);
up = u(:,4:end);
u = u(:,3:end-1);


xi_mm = x_u(1:end-3);
xi_m = x_u(2:end-2);
xi_p = x_u(4:end);
xi = x_u(3:end-1);
xi_u = xi_m;

h_im = xi_m-xi_mm;
h_i  = xi -xi_m;
h_ip = xi_p - xi;
% umm = circshift(u,[0,2]);
% um  = circshift(u,[0 ,1]);
% up  = circshift(u,[0 -1]);
% upp = circshift(u,[0,-2]);

% h_i = circshift(x,-1) - circshift(x,0);
% h_i(end) =L_domain-x(end);%h_i(1);
% h_im  = circshift(h_i,1);
% h_imm = circshift(h_i,2);
% h_ip = circshift(h_i,-1);
% 
% % h_imm = circshift(x,1) - circshift(x,2);
% % h_imm(2)=L_domain-x(end);
% % h_im  = circshift(x,0) - circshift(x,1);
% % h_im(1) = L_domain-x(end);
% % 
% % h_i = circshift(x,-1) - circshift(x,0);
% % h_i(end) =L_domain-x(end);%h_i(1);
% 
% h_ip = circshift(h_i,-1);
% xi_mm = x-h_im;%x - dist_x_pl;
% xi_m  =  x ;
% xi    =  x +h_i;%+ dist_x_pl;
% xi_p  =  x+h_i+h_ip;%x + 2*dist_x_pl;




p0p = umm.*((xiq -xi_m).*(xiq- xi))./((h_im).*(h_im+h_i) )- ...
      um .* ((xiq -xi_mm).*(xiq- xi))./((h_im).*(h_i) )+ ...
      u .* ((xiq -xi_mm).*(xiq- xi_m))./((h_im+h_i).*(h_i) );

p1p = um.*((xiq - xi).*(xiq -xi_p))./((h_i).*(h_i+h_ip))- ...
      u .*((xiq - xi_m).*(xiq - xi_p))./((h_i).*(h_ip))+ ...
      up.*((xiq - xi_m).*(xiq - xi))./((h_i+h_ip).*(h_ip));
%   
% p0p_x = umm .* ((xiq - xi) + (xiq - xi_p))./(dist_x_pl.*(2*dist_x_pl)) - ...
%       um .* ((xiq - xi_m) + (xiq - xi_p))./((dist_x_pl).*(dist_x_pl)) +...
%       u .*( (xiq - xi_m) + (xiq- xi))./((dist_x_pl).*(2* dist_x_pl));
% p1p_x = um.* ((xiq - xi_p) + (xiq - xi_pp))./((dist_x_pl).*(2*dist_x_pl))-...
%       u .* ((xiq - xi) + (xiq - xi_pp))./((dist_x_pl).*(dist_x_pl))+...
%       up .*( (xiq - xi) + (xiq- xi_p))./((2*dist_x_pl).*(dist_x_pl));
% 

p0p_x = umm.*((xiq -xi_m)+(xiq- xi))./((h_im).*(h_im+h_i) )- ...
      um .* ((xiq -xi_mm)+(xiq- xi))./((h_im).*(h_i) )+ ...
      u .* ((xiq -xi_mm)+(xiq- xi_m))./((h_im+h_i).*(h_i) );

p1p_x = um.*((xiq - xi)+(xiq -xi_p))./((h_i).*(h_i+h_ip))- ...
      u .*((xiq - xi_m)+(xiq - xi_p))./((h_i).*(h_ip))+ ...
      up.*((xiq - xi_m)+(xiq - xi))./((h_i+h_ip).*(h_ip));


% p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );

% Smooth Indicators (Beta factors)
dx = x(2)-x(1);
mat = (1/(6*dx)) * [-11, 18, -9, 2;...
                -2, -3 ,6 , -1;...
                 1, -6, 3, 2 ;...
                 -2, 9 , -18, 11];

             
u_arrs_h = [umm(1,:);um(1,:);u(1,:);up(1,:)];
u_ders_h = mat*u_arrs_h;

u_arrs_hv = [umm(2,:);um(2,:);u(2,:);up(2,:)];
u_ders_hv = mat*u_arrs_hv;

u_dmm  = [u_ders_h(1,:);u_ders_hv(1,:)];
u_dm   = [u_ders_h(2,:);u_ders_hv(2,:)];
u_d   = [u_ders_h(3,:);u_ders_hv(3,:)];
u_dp  = [u_ders_h(4,:);u_ders_hv(4,:)];


% B0p  = .25 * (abs(u - umm) - abs(4*um - 3*umm +u)).^2;
% B1p = .25*(abs(u - umm) - abs(4*um - umm -3*u)).^2;

B0p = 4*(abs(u_d - u_dm) - abs(u_dm-u_dmm)).^2;
B1p = 4*(abs(u_dp - u_d) - abs(u_d-u_dm)).^2;
% [nonu_ders, nonu_inds] =  fn_nonu_ders_nonu_non_periodic(x_u,u_ders,L_domain);
% B0p  = nonu_inds(1,:);
% B1p = nonu_inds(2,:);

% linear weights
epsilon = 1e-6;
d0p = - (xiq - xi_p)./(xi_p - xi_mm);

d1p = (xiq - xi_mm)./(xi_p - xi_m);
d0p_x = -1./(xi_p - xi_mm);
d1p_x = 1./(xi_p - xi_mm);
% Alpha weights 
% pow = 1;
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
% poly_rec = d0n.*p0n + d1n.*p1n;
poly_rec_p_x = w0p.*p0p_x + w1p.*p1p_x + w0p_x.*p0p + w1p_x.*p1p;% poly_rec_x = (poly_rec_n-circshift(poly_rec_n,[0,1])+poly_rec_p-circshift(poly_rec_p,[0,1]))/dx;
poly_rec =poly_rec_n + poly_rec_p;
poly_rec_x =( poly_rec_n_x + poly_rec_p_x);% (poly_rec_n-circshift(poly_rec_n,[0,1])+poly_rec_p-circshift(poly_rec_p,[0,1]))/dx;%;



dx_v = xi_v(2:end)-xi_v(1:end-1);
dx_u = xi_u(2:end)-xi_u(1:end-1);
dx_v = [dx_v, dx_v(end)];
dx_u  =[dx_u, dx_u(end)];
% poly_rec_x =(1./dx_v).*[poly_rec_n(:,2:end)-poly_rec_n(:,1:end-1)]+ (1./dx_u).*[poly_rec_p(:,2:end)-poly_rec_p(:,1:end-1),[0;0]];%;

end
