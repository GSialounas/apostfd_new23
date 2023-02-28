function [poly_rec, poly_rec_x] = fn_WENO3_rec_nonu_non_periodic(xiq, x,dist_x_pl, w,L_domain)


% 

%%%%%%+++++++++++++++
% %% flux splitting
 c=1; flux = @(w) c*w;
        dflux = @(w) c*ones(size(w));

a = max(abs(dflux(w)));  v=0.5*(flux(w)+a*w); u=circshift(0.5*(flux(w)-a*w),[0,-1]);
% We are seeking an approximation p(x_{iq}), of the function v(x_{iq}),
% where x_iq is in the interval [x_{i-1}, x_{i}].  We will look at two 
% candidate stencils: 
% {x_{i-2}, x_{i-1} x_{i}} and {x_{i-1}, x_{i}, x_{i+1}}

vmm = circshift(v,[0,2]);
vm  = circshift(v,[0 1]);
vp  = circshift(v,[0 -1]);
vpp = circshift(v,[0,-2]);


h_i = circshift(x,-1) - x;
h_i(end) = L_domain - x(end);

h_imm = circshift(h_i,2);
h_im  =circshift(h_i,1);
h_ip = circshift(h_i,-1);

% h_imm = circshift(x,1) - circshift(x,2);
% h_imm(2)=L_domain-x(end);
% h_im  = circshift(x,0) - circshift(x,1);
% h_im(1) = L_domain-x(end);
% 
% h_i = circshift(x,-1) - circshift(x,0);
% h_i(end) =L_domain-x(end);%h_i(1);

xi_mm = x-h_im-h_imm;%x - 2*dist_x_pl;
xi_m  =  x-h_im;%x - dist_x_pl;
xi    =  x;  
xi_p  =  x+h_i;%x + dist_x_pl;

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

             
v_arrs = [vmm;vm;v;vp];
v_ders = mat*v_arrs;

v_dmm = v_ders(1,:);
v_dm  =v_ders(2,:);
v_d  =v_ders(3,:);
v_dp  =v_ders(4,:);

% B0n = 4*(abs(v_d - v_dm) - abs(v_dm-v_dmm)).^2;
% B1n = 4*(abs(v_dp - v_d) - abs(v_d-v_dm)).^2;

[nonu_ders, nonu_inds] =  fn_nonu_ders(x,v,L_domain);
B0n  = nonu_inds(1,:);
B1n = nonu_inds(2,:);
% linear weights
epsilon = 1e-6;
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
umm = circshift(u,[0,2]);
um  = circshift(u,[0 ,1]);
up  = circshift(u,[0 -1]);
upp = circshift(u,[0,-2]);


h_i = circshift(x,-1) - circshift(x,0);
h_i(end) =L_domain-x(end);%h_i(1);
h_im = circshift(h_i,1);
h_ip = circshift(h_i,-1);
% 
% h_imm = circshift(x,1) - circshift(x,2);
% h_imm(2)=L_domain-x(end);
% h_im  = circshift(x,0) - circshift(x,1);
% h_im(1) = L_domain-x(end);
% 
% h_i = circshift(x,-1) - circshift(x,0);
% h_i(end) =L_domain-x(end);%h_i(1);

h_ip = circshift(h_i,-1);
xi_mm = x-h_im;%x - dist_x_pl;
xi_m  =  x ;
xi    =  x +h_i;%+ dist_x_pl;
xi_p  =  x+h_i+h_ip;%x + 2*dist_x_pl;




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

             
u_arrs = [umm;um;u;up];
u_ders = mat*u_arrs;

u_dmm  = u_ders(1,:);
u_dm   = u_ders(2,:);
u_d   = u_ders(3,:);
u_dp  = u_ders(4,:);

% B0p = 4*(abs(u_d - u_dm) - abs(u_dm-u_dmm)).^2;
% B1p = 4*(abs(u_dp - u_d) - abs(u_d-u_dm)).^2;
[nonu_ders, nonu_inds] =  fn_nonu_ders(x,u,L_domain);
B0p  = nonu_inds(1,:);
B1p = nonu_inds(2,:);

% linear weights
epsilon = 1E-6;
d0p = - (xiq - xi_p)./(xi_p - xi_mm);

d1p = (xiq - xi_mm)./(xi_p - xi_m);
d0p_x = -1./(xi_p - xi_mm);
d1p_x = 1./(xi_p - xi_mm);
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
% poly_rec = d0n.*p0n + d1n.*p1n;
poly_rec_p_x = w0p.*p0p_x + w1p.*p1p_x + w0p_x.*p0p + w1p_x.*p1p;% poly_rec_x = (poly_rec_n-circshift(poly_rec_n,[0,1])+poly_rec_p-circshift(poly_rec_p,[0,1]))/dx;
poly_rec = poly_rec_p +poly_rec_n;
poly_rec_x = (poly_rec_n_x+poly_rec_p_x);% (poly_rec_n-circshift(poly_rec_n,[0,1])+poly_rec_p-circshift(poly_rec_p,[0,1]))/dx;%;


%%%%%-----------

% 
%  c=1; flux = @(w) c*w;
%         dflux = @(w) c*ones(size(w));
% 
% a = max(abs(dflux(w)));  v=0.5*(flux(w)+a*w); u=circshift(0.5*(flux(w)-a*w),[0,-1]);
% % We are seeking an approximation p(x_{iq}), of the function v(x_{iq}),
% % where x_iq is in the interval [x_{i-1}, x_{i}].  We will look at two 
% % candidate stencils: 
% % {x_{i-2}, x_{i-1} x_{i}} and {x_{i-1}, x_{i}, x_{i+1}}
% 
% vmm = circshift(v,[0,2]);
% vm  = circshift(v,[0 1]);
% vp  = circshift(v,[0 -1]);
% vpp = circshift(v,[0,-2]);
% 
% 
% h_i = circshift(x,-1) - x;
% h_i(end) = L_domain - x(end);
% 
% 
% h_im  =circshift(h_i,1);
% h_ip = circshift(h_i,-1);
% h_ipp = circshift(h_i,-2);
% 
% % h_imm = circshift(x,1) - circshift(x,2);
% % h_imm(2)=L_domain-x(end);
% % h_im  = circshift(x,0) - circshift(x,1);
% % h_im(1) = L_domain-x(end);
% % 
% % h_i = circshift(x,-1) - circshift(x,0);
% % h_i(end) =L_domain-x(end);%h_i(1);
% 
% xi_m = x-h_im;%x - 2*dist_x_pl;
% xi    =  x;  
% xi_p  =  x+h_i;%x + dist_x_pl;
% xi_pp = x + h_i+ h_ip;
% 
% % Polynomial with stencil {x_{i-1}, x_i, x_{i+1}}
% 
% p0n = vm.*((xiq - xi).*(xiq - xi_p))./(h_im.*(h_im+h_i)) - ...
%       v.*((xiq-xi_m).*(xiq-xi_p))./(h_im.*h_i)+...
%       vp.*((xiq - xi_m).*(xiq - xi))./((h_im + h_i).*(h_i));
%   
% p1n = v.*((xiq - xi_p).*(xiq - xi_pp))./(h_i.*(h_i+h_ip)) - ...
%       vp.*((xiq-xi).*(xiq-xi_pp))./(h_i.*h_ip)+...
%       vpp.*((xiq - xi).*(xiq - xi_p))./((h_i + h_ip).*(h_ip));
% 
% p0n_x = vm.*((xiq - xi).*(xiq - xi_p))./(h_im.*(h_im+h_i)) - ...
%         v.*((xiq-xi_m).*(xiq-xi_p))./(h_im.*h_i)+...
%         vp.*((xiq - xi_m).*(xiq - xi))./((h_im + h_i).*(h_i));
%   
% p1n_x = v.*((xiq - xi_p) + (xiq - xi_pp))./(h_i.*(h_i+h_ip)) - ...
%         vp.*((xiq-xi) + (xiq-xi_pp))./(h_i.*h_ip)+...
%         vpp.*((xiq - xi) + (xiq - xi_p))./((h_i + h_ip).*(h_ip));
% 
% 
% % p0n = dist_x_pl .* v .* ((xiq - (x-dist_x_pl)) + (xiq -(x+dist_x_pl))) ./ ( (dist_x_pl).*(-dist_x_pl) );
% % p1n = dist_x_pl .* (v + vp) .* ((xiq - (x-dist_x_pl)) + (xiq -x)) ./ ( (2*dist_x_pl).*(dist_x_pl) );
% 
% % Smooth Indicators (Beta factors)
% % dx = x(2)-x(1);
% % mat = (1/(6*dx)) * [-11, 18, -9, 2;...
% %                 -2, -3 ,6 , -1;...
% %                  1, -6, 3, 2 ;...
% %                  -2, 9 , -18, 11];
% % 
% %              
% % v_arrs = [vmm;vm;v;vp];
% % v_ders = mat*v_arrs;
% % 
% % v_dmm = v_ders(1,:);
% % v_dm  =v_ders(2,:);
% % v_d  =v_ders(3,:);
% % v_dp  =v_ders(4,:);
% 
% % B0n = 4*(abs(v_d - v_dm) - abs(v_dm-v_dmm)).^2;
% % B1n = 4*(abs(v_dp - v_d) - abs(v_d-v_dm)).^2;
% 
% [nonu_ders, nonu_inds] =  fn_nonu_ders_nonu(x,v,L_domain);
% B0n  = nonu_inds(1,:);
% B1n = nonu_inds(2,:);
% % linear weights
% epsilon = 1E-6;
% d0n = - (xiq - xi_pp)./(xi_pp - xi_m);
% 
% d1n = (xiq - xi_m)./(xi_pp - xi_m);
% d0n_x = -1./(xi_pp - xi_m);
% d1n_x = 1./(xi_pp - xi_m);
% 
% % Alpha weights 
% pow = 1;
% alpha0n = d0n./(epsilon + B0n).^pow;
% alpha1n = d1n./(epsilon + B1n).^pow;
% alphasumn = alpha0n + alpha1n;
% 
% alpha0n_x = d0n_x./(epsilon + B0n).^pow;
% alpha1n_x = d1n_x./(epsilon + B1n).^pow;
% alphasumn_x =  alpha0n_x +alpha1n_x;
% 
% 
% % ENO stencils weights
% w0n = alpha0n./alphasumn;
% w1n = alpha1n./alphasumn;
% w0n_x = ((alpha0n_x).*(alphasumn) - alpha0n.*alphasumn_x )./ (alphasumn.^2);
% w1n_x = ((alpha1n_x).*(alphasumn) - alpha1n.*alphasumn_x )./ (alphasumn.^2);
% 
% % Numerical value at 
% poly_rec = w0n.*p0n + w1n.*p1n;
% poly_rec_x = w0n.*p0n_x + w1n.*p1n_x + w0n_x.*p0n + w1n_x.*p1n;
