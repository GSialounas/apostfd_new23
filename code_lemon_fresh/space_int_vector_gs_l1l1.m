function int = space_int_vector_gs_l1l1(x,dist_x_pl,dist_x_min,uold, u,tj,ex,c_0_coeff_arr,mask) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq
% c_coefficients = [c_0; c_1; c_2; c_3];
% h = x(2)-x(1);
int = 0;
% for j = 1 : length(x)-1
l_coeffs = length(c_0_coeff_arr);

% % diff_x_arr = zeros(size(c_0_coeff_arr));
% mask = ones(size(x));
% mask(abs(x)>1)=0;
h = x(2)-x(1);
dt= .1*h;
for iqt = 1:nq
    tq(iqt) = 0.5*dt*xq(iqt) + dt/2 +tj;
    diff_t = tq(iqt) - tj;
    c_0_t = uold;
    c_1_t = (1/dt)*(u-uold);
    
    Ut = c_0_t + c_1_t*diff_t;
    c_0_ts = Ut; % this is u at t_n
    c_1_ts = 1./dist_x_pl.*(circshift(Ut,-1)-circshift(Ut,0));%Ut_x; % - f_h(U^n)
%     IU_temp = (tq(iqt) - tj)*u/dt + (dt+tj-tq(iqt))*uold/dt;
    int_x=  0;
    
    for iq = 1:nq
        xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
        diff_x = xiq-x;
        
        IU = c_0_ts + c_1_ts .* diff_x ;
        
        
        %     [IU, IUx] = fn_WENO3_rec_burger(xiq, x,dist_x_pl, u);
        
        
        int_x = int_x + sum(wq(iq).* dist_x_pl .* mask.*abs(IU - ex(xiq,tq(iqt))));
    end
    int  = int+ dt*wq(iqt)*int_x;
end
% end

end