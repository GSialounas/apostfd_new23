function [L2Rt,L2Rt_arr,L1Rt_l1,L1Rt_arr_l1,c_0_coeff_arr_new, c_0_coeff_arr_old,max_IUx,xiq,IU,IUt,IUx,mask] = compute_burger_rec_hermite_time_weno_space_non_periodic(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,f_h_old,f_h_new,spatial_disc,flux_fn,c,u_vec,ex)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element

fluxfun='burgers'; % select flux function
% Define our Flux function
switch fluxfun
    case 'linear'   % Scalar Advection, CFL_max: 0.65
        c=1; flux = @(w) c*w; 
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.40  
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
        c=1; flux_lin = @(w) c*w; 
        dflux_lin = @(w) c*ones(size(w));
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


L2Rt = 0;
L1Rt_l1 = 0;

L2Rt_t = 0;
L1Rt_l1_t = 0;

L2Rt_x = 0;
L1Rt_l1_x = 0;
dx =x(2)- x(1);
%
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients
c_0_t = uold; % this is u at t_n
c_1_t = -f_h_old; % - f_h(U^n)
c_2_t = 3/(dt^2) * (u - uold) + (1/dt) * (f_h_new + 2 * f_h_old);
c_3_t = -1/(dt^2) .* (f_h_new + f_h_old) - 2/(dt^3) .* (u - uold);

diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t *diff_t +c_2_t *diff_t^2 + c_3_t*diff_t^3;
Ut_t = c_1_t + 2*c_2_t*diff_t +3*c_3_t*diff_t.^2;

% We calculate the spatial derivative discretisation using the spatial
% reconstruction at time t
h=dist_x_pl(1);
% if spatial_disc == "CS"
%     Ut_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut,1) + circshift(Ut,-1));
%     Ut_xt = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut_t,1) + circshift(Ut_t,-1));
% elseif spatial_disc == "BS"
%     Ut_x = 1./(dist_x_min).*(-circshift(Ut,1) + circshift(Ut,0));
%     Ut_xt = 1./(dist_x_min).*(-circshift(Ut_t,1) + circshift(Ut_t,0));
% elseif spatial_disc == "2S"
%     Ut_x= 1./(2*h)*( 3*circshift(Ut,0)  - 4*circshift(Ut,1) + circshift(Ut,2) );
%     Ut_xt= 1./(2*h)*( 3*circshift(Ut_t,0)  - 4*circshift(Ut_t,1) + circshift(Ut_t,2) );
% elseif spatial_disc == "3S"
%     Ut_x = 1./(6*h)*( 2*circshift(Ut,-1)  + 3*circshift(Ut,0) - 6* circshift(Ut,1) +circshift(Ut,2)); 
%     Ut_xt = 1./(6*h)*( 2*circshift(Ut_t,-1)  + 3*circshift(Ut_t,0) - 6* circshift(Ut_t,1) +circshift(Ut_t,2));
% elseif spatial_disc == "4CS"
%     Ut_x = -1./(12*h).*(-circshift(Ut,2)+8*circshift(Ut,1) - 8*circshift(Ut,-1) + circshift(Ut,-2));    
%     Ut_xt = -1./(12*h).*(-circshift(Ut_t,2)+8*circshift(Ut_t,1) - 8*circshift(Ut_t,-1) + circshift(Ut_t,-2));
% elseif spatial_disc == "WENO5"
%     Ut_x    =   WENO5resAdv1d_fdm_gs(Ut,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);%
%     Ut_xt   =   WENO5resAdv1d_fdm_gs(Ut_t,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t); %
%     uold_x  =   WENO5resAdv1d_fdm_gs(uold,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);%
%     u_x     =   WENO5resAdv1d_fdm_gs(u,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t); %
% 
% 
% elseif spatial_disc == "WENO3"
%     Ut_x    =   WENO3resAdv1d(Ut,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);%
%     Ut_xt   =   WENO3resAdv1d(Ut_t,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t); %
%     uold_x  =   WENO3resAdv1d(uold,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);%
%     u_x     =   WENO3resAdv1d(u,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t); %
% elseif spatial_disc =="LxF"
%     Ut_x    =   WENO3resAdv1d(Ut,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);%
%     Ut_xt   =   WENO3resAdv1d(Ut_t,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t); %
%     uold_x  =   WENO3resAdv1d(uold,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);%
%     u_x     =   WENO3resAdv1d(u,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t); %
% 
% else
% end

% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
% c_0_ts = Ut; % this is u at t_n
% c_1_ts = Ut_x; % - f_h(U^n)
% c_2_ts = 3./(dist_x_pl.^2) .* (circshift(Ut,-1) - Ut) - (1./dist_x_pl) .* (circshift(Ut_x,-1) + 2 * Ut_x);
% c_3_ts = 1./(dist_x_pl.^2) .* (circshift(Ut_x,-1) + Ut_x) - 2./(dist_x_pl.^3) .* (circshift(Ut,-1) - Ut);
% 
% c_0_ts_t = Ut_t; % this is u at t_n
% c_1_ts_t = Ut_xt; % - f_h(U^n)
% c_2_ts_t = 3./(dist_x_pl.^2) .* (circshift(Ut_t,-1) - Ut_t) - (1./dist_x_pl) .* (circshift(Ut_xt,-1) + 2 * Ut_xt);
% c_3_ts_t = 1./(dist_x_pl.^2) .* (circshift(Ut_xt,-1) + Ut_xt) - 2./(dist_x_pl.^3) .* (circshift(Ut_t,-1) - Ut_t);


uold_x = f_h_old;
u_x = f_h_new;
c_0_old = uold;
c_1_old = uold_x;
c_2_old = 3./(dist_x_pl.^2) .* (circshift(uold,-1) - uold) - (1./dist_x_pl).*(circshift(uold_x,-1) + 2 * uold_x);
c_3_old = 1./(dist_x_pl.^2) .* (circshift(uold_x,-1) + uold_x) - 2./(dist_x_pl.^3) .* (circshift(uold,-1) - uold);


c_0_new = u;
c_1_new = u_x;
c_2_new = 3./(dist_x_pl.^2) .* (circshift(u,-1)-u) - (1./dist_x_pl).*(circshift(u_x,-1) + 2 * u_x);
c_3_new = 1./(dist_x_pl.^2) .* (circshift(u_x,-1) + u_x) - (2./dist_x_pl.^3) .* (circshift(u,-1) - u);

c_0_coeff_arr_old= [c_0_old; c_1_old; c_2_old; c_3_old];
c_0_coeff_arr_new=  [c_0_new; c_1_new; c_2_new; c_3_new];


% if tj<1
%     mask = ones(size(x));
% else
%     mask = abs(x-(sqrt(2)*sqrt(1+tj)-1))>6*h;
% end
mask = ones(size(x));
%     xl= zeros(size(nq));
max_IUx = 0;
L2Rt_arr = zeros(size(x));
L1Rt_arr_l1 = zeros(size(x));



for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
%     
%     IU = c_0_ts + c_1_ts .* diff_x + c_2_ts .* diff_x .^2 + c_3_ts .* diff_x.^3;
%     IUx = c_1_ts + 2* c_2_ts .*diff_x + 3* c_3_ts .*diff_x.^2;
%     IUt = c_0_ts_t + c_1_ts_t.*diff_x +  c_2_ts_t.*diff_x.^2 + c_3_ts_t.*diff_x.^3;
%     
    [IU, IUx] = fn_WENO3_rec_burger_non_periodic(xiq, x,dist_x_pl, Ut);
    u_ex = ex(xiq,evalt);

%     [u_ex, u_x_ex, u_t_ex] = ex(xiq,evalt);
%     u_ex  = ex(xiq,evalt);
%     u_t_ex = ones(size(u_ex));
    fU_burger = flux(Ut);
    [~,fUx] = fn_WENO3_rec_burger_non_periodic(xiq, x,dist_x_pl, fU_burger);

    IUt = fn_WENO3_rec_burger_non_periodic(xiq, x,dist_x_pl, Ut_t);

%     R_h_burger= mask.*(IUt + IU.*IUx);
%     L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*mask.*(IUt + IU.*IUx).^2);
%     L2Rt_arr = L2Rt_arr +wq(iq)*dist_x_pl.*mask.*(IUt + IU.*IUx).^2;
%     L1Rt_l1 = L1Rt_l1 + sum(wq(iq)*dist_x_pl.*mask.*abs(IUt + IU.*IUx));
%     L1Rt_arr_l1 = L1Rt_arr_l1 +wq(iq)*dist_x_pl.*mask.*abs(IUt + IU.*IUx);
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*mask.*(IUt + fUx).^2);
    L2Rt_arr = L2Rt_arr +wq(iq)*dist_x_pl.*mask.*(IUt + fUx).^2;
    L1Rt_l1 = L1Rt_l1 + sum(wq(iq)*dist_x_pl.*mask.*abs(IUt + fUx));
    L1Rt_arr_l1 = L1Rt_arr_l1 +wq(iq)*dist_x_pl.*mask.*abs(IUt + fUx);

%     % time derivative error
%     L2Rt_t = L2Rt_t + sum(wq(iq)*dist_x_pl.*mask.*(IUt-u_t_ex).^2);
%     L2Rt_arr_t = L2Rt_arr_t +wq(iq)*dist_x_pl.*mask.*(IUt-u_t_ex).^2;
% 
%     L1Rt_l1_t = L1Rt_l1_t + sum(wq(iq)*dist_x_pl.*mask.*abs(IUt-u_t_ex));
%     L1Rt_arr_l1_t = L1Rt_arr_l1_t +wq(iq)*dist_x_pl.*mask.*abs(IUt-u_t_ex);
%     
%     % space derivative
%     L2Rt_x = L2Rt_x + sum(wq(iq)*dist_x_pl.*mask.*(fUx-fUx_ex).^2);
%     L2Rt_arr_x = L2Rt_arr_x +wq(iq)*dist_x_pl.*mask.*(fUx -fUx_ex).^2;
% 
%     L1Rt_l1_x = L1Rt_l1_x + sum(wq(iq)*dist_x_pl.*mask.*abs(fUx -fUx_ex));
%     L1Rt_arr_l1_x = L1Rt_arr_l1_x +wq(iq)*dist_x_pl.*mask.*abs(fUx -fUx_ex);

    max_IUx = max(max_IUx,max(abs(IUx))); % L_infty norm for IUx at gauss points
end
% end

end
