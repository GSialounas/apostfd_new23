clc;
clear;
close all;

x= linspace(-pi,pi,501);
x=x(1:end-1);
t= 0;
it = 0;
T =2;
fps =10;
showplot=1;

dx = x(2)-x(1);
dt = .1*dx;
uold = zeros(size(x));
% uold(abs(x)<=1)=1;
% uold(x<=-1)=1;
% uold(x>1) = 1;
xl=-1;
xr=1;
ul=1;
ur=1;
ex = @(x,t) fn_exact_step_combo_burger(x,t,ul,ur,xl,xr);
uold= ex(x,0);
u = uold;

figure;



fluxfun='burgers'; % select flux function
% Define our Flux function
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
end

sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        S = @(w) 0.1*w.^2;
    case 'dont'
        S = @(w) zeros(size(w));
end


while t<=T-dt/2
    
    u = uold-dt*.5*(1/dx)*((uold).^2-circshift(uold,1).^2);
    
%     f_h_old  = WENO3_naive_periodic(uold,flux,dflux,S,dx);%resWENO5(uold_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu); %WENO5resAdv1d_fdm_gs(uold_nu,flux,dflux,S,dx);%
%     
%     % Three stage
%     ustage_0 = uold;
%     
%     ustage_1 = uold - dt*f_h_old;
%     f_h_stage_1 =  WENO3_naive_periodic(ustage_1,flux,dflux,S,dx);%resWENO5(ustage_1,c,dx);% fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_1); %WENO5resAdv1d_fdm_gs(ustage_1,flux,dflux,S,dx);%
%     
%     ustage_2 = (.75)*uold +.25*ustage_1 - .25 *dt *(f_h_stage_1);
%     f_h_stage_2 = WENO3_naive_periodic(ustage_2,flux,dflux,S,dx);%resWENO5(ustage_2,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_2); %WENO5resAdv1d_fdm_gs(ustage_2,flux,dflux,S,dx);%
%     
%     ustage_3 = (1/3)*uold +(2/3)*ustage_2 - (2/3)*dt*(f_h_stage_2);
%     
%     u = ustage_3;
%     [f_h_new,f_u] = WENO3_naive_periodic(u,flux,dflux,S,dx);%resWENO5(u_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_3); %WENO5resAdv1d_fdm_gs(u_nu,flux,dflux,S,dx);%
%                 
    t=t+dt;
    it=it+1;
    
    if (showplot &&mod(it-1,fps)==0)
        plot(x,u,'r',x,ex(x,t+dt),'b--')
        axis([-pi pi -.5 1.1])
        title(['t=',num2str(t)])
        legend('u_{approx}','u_{ex}','location','SouthWest')
        xlabel('x')
        ylabel('u')
        pause(0.01)
    end
    uold=u;
end