function [L2Rt,L2Rt_h,L2Rt_hv,L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv,IUt_h, IUt_hv,L2Rt_h_arr,L2Rt_hv_arr, max_IUx] = compute_Rs_vector_temp_1_spatiotemp_1_shw_dambreak(x,dist_x_pl,dist_x_min,h_old,h_new,evalt,tj,dt,f_h_old,f_h_new,spatial_disc,flux_fn)
global nq xq wq

g=9.81;
dx= dist_x_pl(1);
L2Rt = 0;
L2Rt_h = 0;

L2Rt_hv = 0;

L2Rt_arr = zeros(1,length(x)-1);
L2Rt_h_arr = zeros(1,length(x)-1);
L2Rt_hv_arr = zeros(1,length(x)-1);


c_0_coeff_arr_new=0; c_0_coeff_arr_old=0;


c_0_t = h_old;
c_1_t = 1/dt*(h_new-h_old);



diff_t=evalt-tj;
dist_x_pl = x(2:end)-x(1:end-1);
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
ht = c_0_t + c_1_t * diff_t ;
ht_t = c_1_t ;

c_0_ts = ht(:, 1:end-1);
% c_1_ts =  1/dx * (circshift(ht,-1,2)-ht); % periodic
c_1_ts =  [1./dist_x_pl;1./dist_x_pl] .* (ht(:, 2:end) - ht(:, 1:end-1)); % non-periodic



c_0_ts_t = ht_t(:, 1:end-1);
% c_1_ts_t = 1/dx * (circshift(ht_t,-1,2)-ht_t);
c_1_ts_t = [1./dist_x_pl;1./dist_x_pl] .* (ht_t(:, 2:end)-ht_t(:,1:end-1));




%Now we separate the two reconstructions
c_0_ts_h = c_0_ts(1,:);
c_0_ts_hv = c_0_ts(2,:);

c_1_ts_h = c_1_ts(1,:);
c_1_ts_hv = c_1_ts(2,:);


c_0_ts_t_h = c_0_ts_t(1,:);
c_0_ts_t_hv = c_0_ts_t(2,:);

c_1_ts_t_h = c_1_ts_t(1,:);
c_1_ts_t_hv = c_1_ts_t(2,:);



max_IUx_h  = 0;
max_IUx_hv = 0;

for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x(1:end-1) + dist_x_pl/2;
    diff_x = xiq-x(1:end-1);
    
    IU_h   = c_0_ts_h + c_1_ts_h .* diff_x ;
    IU_hv  = c_0_ts_hv + c_1_ts_hv .* diff_x ;
    IUx_h  =  c_1_ts_h ;
    IUx_hv =  c_1_ts_hv ;
    
    IUt_h =  c_0_ts_t_h + c_1_ts_t_h .* diff_x ;
    IUt_hv =  c_0_ts_t_hv + c_1_ts_t_hv .* diff_x ;
    
    R_h = IUt_h +IUx_hv;
    R_hv = IUt_hv + ((IU_h).^(-2)).*(2*IU_hv.*IUx_hv.*IU_h - (IU_hv.^2).*IUx_h) + g*IU_h.*IUx_h;
    
    integral_txq =  (wq(iq)*dist_x_pl.*(R_h).^2) +(wq(iq)*dist_x_pl.*(R_hv).^2);
    integral_txq_h =  (wq(iq)*dist_x_pl.*(R_h).^2);
    integral_txq_hv =  (wq(iq)*dist_x_pl.*(R_hv).^2);

   
    max_IUx_h  =max(max_IUx_h, max(abs(IUx_h)));
    max_IUx_hv = max(max_IUx_hv, max(abs(IUx_hv)));
    
    L2Rt = L2Rt + sum(integral_txq);
    L2Rt_h = L2Rt_h + sum(integral_txq_h);
    L2Rt_hv = L2Rt_hv + sum(integral_txq_hv);
    L2Rt_arr = L2Rt_arr + (integral_txq);
    
    L2Rt_h_arr = L2Rt_h_arr + (integral_txq_h);
    L2Rt_hv_arr = L2Rt_hv_arr + (integral_txq_hv);
    
    
end

max_IUx = max(max_IUx_h, max_IUx_hv);

end