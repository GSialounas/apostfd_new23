function [L2Rt,L2Rt_arr,L2Rt_l1,L2Rt_arr_l1,c_0_coeff_arr_new, c_0_coeff_arr_old,max_IUx,xiq,IU,IUt,IUx,mask] = compute_Rs_vector_temp_1_spatiotemp_1_burger_fresh_non_periodic(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc,flux_fn,ex)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element

% uold = ex(x,tj);
% u = ex(x,tj+dt);

L2Rt = 0;
L2Rt_arr = zeros(size(x));
L2Rt_l1 = 0;
L2Rt_arr_l1 = zeros(size(x));

L2Rt_arr_t = zeros(size(x));
L2Rt_arr_x = zeros(size(x));
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients

% reconstruction with derivative matching at t_n
c_0_t = uold; % this is u at t_n
c_1_t = (1/dt)*(u-uold);



diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t * diff_t ;
Ut_t = c_1_t;



% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
% Reconstruction where we match the spatial derivative at x_j
c_0_ts = Ut(1:end-1); % this is u at t_n
c_1_ts = 1./dist_x_pl.*(Ut(2:end)-Ut(1:end-1));%Ut_x; % - f_h(U^n)

% c_2_ts = (1./(dist_x_pl.^2)) .* (circshift(Ut,-1) - (Ut + dist_x_pl .* Ut_x));

c_0_ts_t = Ut_t(1:end-1); % this is u at t_n
c_1_ts_t = 1./dist_x_pl.*(Ut(2:end)-Ut(1:end-1));%Ut_xt; % - f_h(U^n)

c_0_old = uold(1:end-1);
c_1_old = 1./dist_x_pl.*(uold(2:end)-uold(1:end-1));%uold_x;



c_0_new = u(1:end-1);
c_1_new = 1./dist_x_pl.*(u(2:end)-u(1:end-1));


c_0_coeff_arr_old  =  [c_0_old; c_1_old];
c_0_coeff_arr_new  =  [c_0_new; c_1_new];

h=dist_x_pl(1);
% shock mask
% if tj<1
%     mask = ones(size(x));
% else
%     mask = abs(x-(sqrt(2)*sqrt(1+tj)-1))>20*h;
% end
mask = abs(x-.5*tj)>3*h;
% mask = ones(size(x));

max_IUx = 0;
for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU  = c_0_ts + c_1_ts .* diff_x ;
    IUx =  c_1_ts ;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x;
        FUx = (c_0_ts+c_1_ts.*diff_x).*(c_1_ts);
%     
% [~, fUx] = fn_WENO3_rec_burger(xiq, x,dist_x_pl, .5*Ut.^2);
%     [IU, IUx] = fn_WENO3_rec_burger(xiq, x,dist_x_pl, Ut);
%     IUt = fn_WENO3_rec_burger(xiq, x,dist_x_pl, Ut_t);
    
%     IUt = fn_WENO3_rec_burger(xiq, x,dist_x_pl, Ut_t);
    
    
%     integral_txq = (wq(iq)*dist_x_pl.*mask.*( (IUt+ IU.*IUx)).^2);
%     integral_txq_l1 = (wq(iq)*dist_x_pl.*mask.*( abs(IUt+IU.*IUx)));
    
            integral_txq = (wq(iq)*dist_x_pl.*mask.*( (IUt+FUx)).^2);
    integral_txq_l1 = (wq(iq)*dist_x_pl.*mask.*( abs(IUt+FUx)));
%     integral_txq = (wq(iq)*dist_x_pl.*( (IUt+fUx)).^2);
%     integral_txq_l1 = (wq(iq)*dist_x_pl.*( abs(IUt+fUx)));

    integral_txq_t = (wq(iq)*dist_x_pl.*( IUt).^2);
    integral_txq_x = (wq(iq)*dist_x_pl.*( IUx).^2);
    
    L2Rt = L2Rt + sum(mask.*integral_txq);
    L2Rt_arr = L2Rt_arr + mask.*(integral_txq);
    
    L2Rt_l1 = L2Rt_l1 + sum(mask.*integral_txq_l1);
    L2Rt_arr_l1 = L2Rt_arr_l1 + (mask.*integral_txq_l1);
    
    L2Rt_arr_t =  c_1_ts;%c_2_ts_t.*diff_x.^2 ;% L2Rt_arr_t + (integral_txq_t);
    L2Rt_arr_x =  c_1_ts_t; %c_2_ts_t.*diff_x.^2 ;%L2Rt_arr_x + (integral_txq_x);
    max_IUx = max(max_IUx,max(abs(IUx))); % L_infty norm for IUx at gauss points

    
end
% end

end