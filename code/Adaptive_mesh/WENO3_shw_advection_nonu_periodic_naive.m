function [poly_rec, poly_rec_x] = WENO3_shw_dambreak_nonu_naive(w,flux,dflux,S,dx,eig_Jf_1,eig_Jf_2,F1,F2,x,t,L_domain)

xiq = x;
if length(w)==1
    a=max(abs(dflux(w)));
else
    eig_max_1 = max(abs(eig_Jf_1(w)));
    eig_max_2 = max(abs(eig_Jf_2(w)));
    a= max(eig_max_1,eig_max_2);
end

v=0.5*(flux(w) + a*w);  u=circshift(0.5*(flux(w)-a*w),[0,-1]);

vm  = circshift(v,[0 1]);
vp  = circshift(v,[0 -1]);

h_i = circshift(x,-1)-x;

h_im = circshift(h_i,1);
h_ip = circshift(h_i,-1);

xi = x;
xi_m = xi - h_im;
xi_p = xi+ h_i;
xi_pp = xi+ h_i+h_ip;


p0n = (vm.*h_im).*((xiq-xi_m) + (xiq-xi_p))./(h_im.*(-h_i))+...
      (vm.*h_im +v.*h_i).*((xiq - xi_m) + (xiq - xi))./((h_im + h_i).*(h_i));

p1n =  (v.*h_i).*((xiq-xi) +  (xiq-xi_pp))./(h_i.*(-h_ip))+...
        (v.*h_i + vp.*h_ip) .*((xiq - xi) + (xiq - xi_p))./((h_i + h_ip).*(h_ip));

p0n_x =     (vm.*h_im).*(2)./(h_im.*(-h_i))+...
            (vm.*h_im +v.*h_i).*(2)./((h_im + h_i).*(h_i));
p1n_x = (v.*h_i).*(2)./(h_i.*(-h_ip))+...
        (v.*h_i + vp.*h_ip) .* (2)./((h_i + h_ip).*(h_ip));

epsilon = 1E-6;
d0n = - (xiq - xi_pp)./(xi_pp - xi_m);

d1n = (xiq - xi_m)./(xi_pp - xi_m);
d0n_x = -1./(xi_pp - xi_m);
d1n_x = 1./(xi_pp - xi_m);

B0n = (p0n_x.*h_i).^2;%(vm-v).^2; 
B1n = (p1n_x.*h_i).^2;%(v-vp).^2;
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
poly_rec_n_x = (1./h_im(1:end-1)).*(poly_rec_n(:,2:end)- poly_rec_n(:,1:end-1));
% poly_rec_n_x = (poly_rec_n(2:end)

%%%  

um  = circshift(u,[0 1]);
up  = circshift(u,[0 -1]);

h_i = circshift(h_i,-1); % the h_i on the rhs is the h_i for v
h_im = circshift(h_i,1);
h_ip = circshift(h_i,-1);

xi = [x(2:end),L_domain];
xi_m = xi - h_im;
xi_p = xi+ h_i;
xi_pp = xi+ h_i+h_ip;

p0p = (um.*h_im).*((xiq-xi_m) + (xiq-xi_p))./(h_im.*(-h_i))+...
      (um.*h_im +u.*h_i).*((xiq - xi_m) + (xiq - xi))./((h_im + h_i).*(h_i));

p1p =  (u.*h_i).*((xiq-xi) +  (xiq-xi_pp))./(h_i.*(-h_ip))+...
        (u.*h_i + up.*h_ip) .* ((xiq - xi) + (xiq - xi_p))./((h_i + h_ip).*(h_ip));
    
p0p_x = (um.*h_im).*(2)./(h_im.*(-h_i))+...
        (um.*h_im +u.*h_i).*(2)./((h_im + h_i).*(h_i));
    
p1p_x =  (u.*h_i).*(2)./(h_i.*(-h_ip))+...
        (u.*h_i + up.*h_ip) .* (2)./((h_i + h_ip).*(h_ip));
    
    
    
d0p = - (xiq - xi_pp)./(xi_pp - xi_m);
d1p = (xiq - xi_m)./(xi_pp - xi_m);
d0p_x = -1./(xi_pp - xi_m);
d1p_x = 1./(xi_pp - xi_m);

B0p = (p0p_x.*h_im).^2;%(um-u).^2; 
B1p = (p1p_x.*h_im).^2;%(u-up).^2;

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
poly_rec_p_x = (1./h_im(1:end-1)).*(poly_rec_p(:,2:end)- poly_rec_p(:,1:end-1));
% poly_rec_n_x = (poly_rec_n(2:end)

poly_rec =  poly_rec_n + poly_rec_p;
poly_rec_x  =poly_rec_p_x +poly_rec_n_x;

% poly_rec_x = poly_rec_x(:,1:end-1);
% hn = poly_rec_n;
% hp = poly_rec_p;
% res = (hp(:,2:end)-hp(:,1:end-1)+hn(:,2:end)-hn(:,1:end-1))./h_im;% - S(w);


end