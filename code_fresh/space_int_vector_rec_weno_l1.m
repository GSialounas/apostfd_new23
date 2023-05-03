function [int,int_x,int_t] = space_int_vector_rec_weno_l1(x,dist_x_pl,dist_x_min,u,evalt,ex) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq
 
int = 0;

mask = ones(size(x));
mask(abs(x)>1) = 0;
for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
    diff_x = xiq-x;
%     IU = sum(c_0_coeff_arr.*[ones(1,l_coeffs); diff_x; diff_x.^2; diff_x.^3],1);    
    IU = fn_WENO3_rec(xiq, x,dist_x_pl, u);
    int = int + sum(wq(iq).* mask.* dist_x_pl .* abs(IU - ex(xiq,evalt)));
end

end
