function [L1Rt, L1Rt_arr] = compute_temp_1_space_weno(x,dist_x_pl,uold,u,evalt,tj,dt,mask)
global nq xq wq

c0 = uold ;
c1 = (1/dt) *(u-uold);

diff_t =evalt-tj;

L1Rt = 0;
L1Rt_arr = zeros(size(x));

Ut = c0 + c1*diff_t;
Ut_t = c1;

% mask = ones(size(x));
% mask(abs(x)>1)=0;
h = x(2)-x(1);
for iq = 1:nq
    xiq = .5*dist_x_pl*xq(iq) + x +dist_x_pl/2;
    [IU,IUx] = fn_WENO3_rec(xiq,x,h,Ut);
    [IU_t,~] = fn_WENO3_rec(xiq,x,h,Ut_t);
    [~, FUx] = fn_WENO3_rec(xiq,x,h,.5*Ut.^2);
    
    L1Rt_arr = L1Rt_arr + (wq(iq)*h*mask.*abs(IU_t+FUx));
    L1Rt = L1Rt + sum(wq(iq)*h*mask.*abs(IU_t+FUx));

end

end