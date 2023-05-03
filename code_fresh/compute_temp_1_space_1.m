function [L1Rt, L1Rt_arr, L2Rt, L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old, max_IUx] = compute_temp_1_space_1(x,dist_x_pl,uold,u,evalt,tj,dt,mask)
global nq xq wq

c0 = uold ;
c1 = (1/dt) *(u-uold);

diff_t =evalt-tj;

L1Rt = 0;
L1Rt_arr = zeros(size(x));

L2Rt = 0;
L2Rt_arr = zeros(size(x));

Ut = c0 + c1*diff_t;
Ut_t = c1;

% Spatial reconstruction
c_0_ts = Ut; 
c_1_ts = 1./dist_x_pl.*(circshift(Ut,-1)-circshift(Ut,0));


c_0_ts_t = Ut_t; 
c_1_ts_t = 1./dist_x_pl.*(circshift(Ut_t,-1)-circshift(Ut_t,0));

c_0_old = uold;
c_1_old = 1./dist_x_pl.*(circshift(uold,-1)-uold);%uold_x;


c_0_new = u;
c_1_new = 1./dist_x_pl.*(circshift(u,-1)-u);

c_0_coeff_arr_old  =  [c_0_old; c_1_old];
c_0_coeff_arr_new  =  [c_0_new; c_1_new];


% mask = ones(size(x));
% mask(abs(x)>1)=0;
h = x(2)-x(1);
max_IUx = 0;
for iq = 1:nq
    xiq = .5*dist_x_pl*xq(iq) + x +dist_x_pl/2;
    diff_x = xiq-x;
    
    IU  = c_0_ts + c_1_ts .* diff_x ;
    FUx = (c_0_ts+c_1_ts.*diff_x).*(c_1_ts);
    IU_x =  c_1_ts ;
    IU_t = c_0_ts_t + c_1_ts_t.*diff_x;
    
    L1Rt_arr = L1Rt_arr + (wq(iq)*h*mask.*abs(IU_t+FUx));
    L1Rt = L1Rt + sum(wq(iq)*h*mask.*abs(IU_t+FUx));
    
    L2Rt_arr = L2Rt_arr + (wq(iq)*h*mask.*(IU_t+FUx).^2);
    L2Rt = L2Rt + sum(wq(iq)*h*mask.*(IU_t+FUx).^2);
    max_IUx = max(max_IUx, max(abs(IU_x)));

end

end