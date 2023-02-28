clc;
clear;
close all;


L_domain= 32*pi;
x = linspace(0,L_domain,501);
ratio_cTf = 1;
dx_fine = x(2)-x(1);
part_fine = .5;
x = create_grid(L_domain, ratio_cTf, dx_fine, part_fine);

dist_x_pl = x(2:end)- x(1:end-1);
% dt = .05 *(x(2)-x(1));
dx = x(2)-x(1);
x= x(1:end-1);
fluxfun='linear'; % select flux function
% Define our Flux functionre
switch fluxfun
    case 'linear'   % Scalar Advection, CFL_max: 0.65
        c=1; flux = @(w) c*w;
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.40
        flux = @(w) w.^2/2;
        dflux = @(w) w;
    case 'buckley' % Buckley-Leverett, CFL_max: 0.20 & tEnd: 0.40
        flux = @(w) 4*w.^2./(4*w.^2+(1-w).^2);
        dflux = @(w) 8*w.*(1-w)./(5*w.^2-2*w+1).^2;
    case 'shallow_water'
        % w = [h; hv]
        flux = @(w) [w(2,:); (w(2,:).^2./w(1,:))+.5*g*w(1,:).^2];
        eig_Jf_1 = @(w) w(2,:)./w(1,:) + sqrt(g*w(1,:));
        eig_Jf_2 = @(w) w(2,:)./w(1,:) - sqrt(g*w(1,:));
        dflux = @(w) [0 , 1; -(w(2,:).^2./(w(1,:).^2)) + g*w(1,:), 2*w(2)./w(1)];
end


sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        %         S = @(w) 0.1*w.^2;
        
        %          S = @(w) .1*w.^2;
        %         S =@(x,t) fn_exact_forcing(x,t,d,a,omega,k,g);
    case 'dont'
        S = @(w) zeros(size(w));
end
% x= x(1:end-1);
w = zeros(size(x));
w(x>10 & x<20) = 1;

% w=sin(2*pi*x);
% figure
% plot(x,w);
dt= .1*dx;
t=0;
T = 20000*dt;
w_new = zeros(size(w));
it=0;
figure;
while t<T
    
    [w_rec, f_h] =  fn_WENO3_rec_shw_trial_n2(x, x,dist_x_pl, w,flux,dflux,L_domain);%fn_WENO3_rec_shw_nonu(x, x,dist_x_pl, w,flux,dflux);%fn_WENO3_rec_shw(x, x,dist_x_pl, w);%
    f_h =WENO3resAdv1d(w,flux,dflux,S,dx);

%     f_h = WENO3resAdv1d(w,flux,dflux,S,dx);%fn_WENO3_rec_shw_trial(x, x,dist_x_pl, w,flux,dflux);%WENO3resAdv1d_non_periodic(w,flux,dflux,S,dx);
    w_new = w -dt * f_h;%WENO3resAdv1d_non_periodic(w,flux,dflux,S,dx);
    
    w = w_new;
    
    if mod(it-1,20)==0
%         plot(x,w,'r.')

        plot(x,w,'b',x,w_rec,'r.')
        grid on;
        pause(0.01)
    end
    t =t +dt;
    it = it+1;
end

% x = 1:10;
% v= ones(size(x));
% dist_x_pl  = ones(size(x));
% vm = circshift(v,1);
% vmm = circshift(v,2);
% vp = circshift(v,-1);
% vpp = circshift(v,-2);
% 
% dx = x(2) - x(1);
% xi = x;
% xi_mm = x- dx;
% xi_m = x - dx;
% xi_p = x + dx;
% xi_pp = x  + 2*dx;
% xiq = xi_pp;


% 
% pol = vm.* (xiq - xi).*(xiq - xi_p).*(xiq - xi_pp)./((-dist_x_pl).*(-2*dist_x_pl).*(-3*dist_x_pl))+...
%       v .* (xiq - xi_m).*(xiq - xi_p).*(xiq - xi_pp)./((dist_x_pl).*(-dist_x_pl).*(-2*dist_x_pl))+...
%       vp .* (xiq - xi_m).*(xiq- xi).*(xiq- xi_pp)./((2*dist_x_pl).*(dist_x_pl).*(-dist_x_pl))+...
%       vpp .* (xiq - xi_m).*(xiq- xi).*(xiq- xi_p)./((3*dist_x_pl).*(2*dist_x_pl).*(dist_x_pl));
% pol_x = vm.*( (xiq - xi_p).*(xiq - xi_pp) +  (xiq - xi).*(xiq - xi_pp)+  (xiq - xi).*(xiq - xi_p))./((-dist_x_pl).*(-2*dist_x_pl).*(-3*dist_x_pl))+...
%         v .* ((xiq - xi_p).*(xiq - xi_pp) +  (xiq - xi_m).*(xiq - xi_pp)  +  (xiq - xi_m).*(xiq - xi_p))./((dist_x_pl).*(-dist_x_pl).*(-2*dist_x_pl))+...
%         vp .* ((xiq- xi).*(xiq- xi_pp) +  (xiq - xi_m).*(xiq- xi_pp)  +  (xiq - xi_m).*(xiq- xi))./((2*dist_x_pl).*(dist_x_pl).*(-dist_x_pl))+...
%         vpp .* ((xiq- xi).*(xiq- xi_p) + (xiq - xi_m).*(xiq- xi_p)  +  (xiq - xi_m).*(xiq- xi) )./((3*dist_x_pl).*(2*dist_x_pl).*(dist_x_pl));
% 
% 
% 
% pol_2 = vmm .* (xiq - xi_m).*(xiq - xi).*(xiq - xi_p)./((-dist_x_pl).*(-2*dist_x_pl).*(-3*dist_x_pl))+...
%         vm  .* (xiq - xi_mm).*(xiq - xi).*(xiq - xi_p)./((dist_x_pl).*(-dist_x_pl).*(-2*dist_x_pl))+...
%         v   .* (xiq - xi_mm).*(xiq- xi_m).*(xiq- xi_p)./((2*dist_x_pl).*(dist_x_pl).*(-dist_x_pl))+...
%         vp  .* (xiq - xi_mm).*(xiq- xi_m).*(xiq- xi)./((3*dist_x_pl).*(2*dist_x_pl).*(dist_x_pl));
%     
% pol_2_x = vmm .*(xiq - xi_m).*(xiq - xi).*(xiq - xi_p)./((-dist_x_pl).*(-2*dist_x_pl).*(-3*dist_x_pl))+...
%         vm  .* (xiq - xi_mm).*(xiq - xi).*(xiq - xi_p)./((dist_x_pl).*(-dist_x_pl).*(-2*dist_x_pl))+...
%         v   .* (xiq - xi_mm).*(xiq- xi_m).*(xiq- xi_p)./((2*dist_x_pl).*(dist_x_pl).*(-dist_x_pl))+...
%         vp  .* (xiq - xi_mm).*(xiq- xi_m).*(xiq- xi)./((3*dist_x_pl).*(2*dist_x_pl).*(dist_x_pl));
% 
% 
% 

