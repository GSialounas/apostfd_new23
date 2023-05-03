function [int_l1, int_l2] = space_int_vector_gs_l1_non_periodic(x,dist_x_pl,dist_x_min,u,evalt,ex,c_0_coeff_arr) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq
% c_coefficients = [c_0; c_1; c_2; c_3];
% h = x(2)-x(1);
int_l1 = 0;
int_l2 = 0;

% for j = 1 : length(x)-1
l_coeffs = length(c_0_coeff_arr);

% diff_x_arr = zeros(size(c_0_coeff_arr));
for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
    diff_x = xiq-x;
    
%     [IU, IUx] = fn_WENO3_rec_burger(xiq, x,dist_x_pl, u);
    [IU, IUx]  = fn_WENO3_rec_burger_non_periodic(xiq, x,dist_x_pl, u);
%     IU = sum(c_0_coeff_arr.*(diff_x.^([0:(size(c_0_coeff_arr,1)-1)]')),1);
    
    int_l1 = int_l1 + sum(wq(iq).* dist_x_pl .* abs(IU - ex(xiq,evalt)));
    int_l2 = int_l2 + sum(wq(iq).* dist_x_pl .* abs(IU - ex(xiq,evalt)));
end
% end

end