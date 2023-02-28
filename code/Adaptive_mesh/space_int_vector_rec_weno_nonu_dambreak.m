function int = space_int_vector_rec_weno_nonu_dambreak(x,dist_x_pl,dist_x_min,u,evalt,ex,L_domain) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq
 
int = 0;
for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x(1:end-1) + .5*dist_x_pl;
    diff_x = xiq-x(1:end-1);
%     IU = sum(c_0_coeff_arr.*[ones(1,l_coeffs); diff_x; diff_x.^2; diff_x.^3],1);    
%     IU = fn_WENO3_rec_nonu(xiq, x,dist_x_pl, u,L_domain);
%         IU = fn_WENO3_rec_nonu_onefrompaper(xiq, x,dist_x_pl, u,L_domain);
    [IU, IUx] = fn_WENO3_rec_nonu_dambreak_error_only(xiq, x(1:end-1),dist_x_pl, [u(1,:);u(2,:)],L_domain);
    [ex_h, ex_hv] = ex(xiq,evalt);
%     int = int + sum(wq(iq).* dist_x_pl .* (IU - ex(xiq,evalt)).^2);
%     int = int + sum(wq(iq).* dist_x_pl .* (IU - ex(xiq,evalt)).^2);
    int = int + sum(wq(iq).* dist_x_pl .* ((IU(1,:) - ex_h).^2  + (IU(2,:) - ex_hv).^2));
end

end
