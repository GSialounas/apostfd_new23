clc;
clear;
close all;

fps = 50;
showplot = 1;
global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();
scheme_str = 'FTBS';
init_cond = 'hatIC';

ul = 1;
ur = 1;

if init_cond == 'stpIC'
    ex_step = @(x,t) fn_exact_step_burger(x,t,ul,ur);
    ex=ex_step;
    ex_step_x = @(x,t) zeros(size(x));
elseif init_cond =='hatIC'
    ex = @(x,t) fn_burger_ex_lin(x,t);
elseif init_cond == 'cmbIC'
    xl=-1;
    xr=1;
    ex = @(x,t) fn_exact_step_combo_burger(x,t,ul,ur,xl,xr);
else
    ex_smooth  = @(x,t) sin(10*pi*x);
    ex_smooth_x = @(x,t) 10*pi*cos(10*pi*x);
    ex= ex_smooth;
end


cell_arr = {};


for m = 1:4
    disp(m)
    h = 2^(-m-8);
     x=0:h:1;
%     x = -pi:h:pi;
    x= -pi+2*pi*(x-x(1))/(x(end)-x(1));
    len_x = length(x);
%     x = linspace(-2*pi,2*pi,len_x+3);
    dist_x_pl =  x(2:end) -x(1:end-1);
    dist_x_min = dist_x_pl;
    h = x(2)-x(1);
    dt= .1*h;
    if m ==1
        if init_cond == 'stpIC'
            T_multiples = 1;
            N_steps = 800;
            T = T_multiples*N_steps*dt;
        elseif init_cond == 'hatIC'
            T_multiples = 3;
            N_steps =800;
            T= T_multiples*N_steps*dt;
        elseif init_cond == 'cmbIC'
            T_multiples= 1;
            N_steps = 800;
            T =  T_multiples*N_steps*dt;
        end
    end
   
    it =0;
    t = 0;

    % ic


    
    L1L1R = 0;
    L2L2R = 0;
    exponential_factor = 0;

    h_arr(m) = h;
    x=x(1:end-1);
    
    mask = ones(size(x));
    if init_cond == 'stpIC'
        mask(abs(x)>1) =0;
    elseif init_cond == 'hatIC'
    elseif init_cond == 'cmbIC'
        mask(abs(x)>2.0) = 0;
    else
    end

    uold = ex(x,t);
    u = uold;
    bound_arr = [];
    error_arr = [];
    bound_arr_l2 = [];
    error_arr_l2 = [];
    error_arr_l1l1 = []; %L1 time L1 space
    EI_arr = [];
    EI_arr_ohl = [];

    time_arr= [0];
    etat_arr = [];
    etac_arr = [];
    bound_ohl_arr = [];
    etat = 0;
    etac = 0;
    bound_ohl = 0;
    error_0 = 0;
    
%     for iq = 1:nq
%         xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
%         diff_x = xiq-x;
%         IU = fn_WENO3_rec(xiq, x,dist_x_pl, uold);        
%         
%         error_0 = error_0 + sum(wq(iq).* dist_x_pl .* abs(IU - ex_step(xiq,0)));
%     end
    % this is error 0  using linear rec
%     [error_0, error_0_l2] = space_int_vector_gs_l1(x,dist_x_pl, dist_x_min,u, 0,ex,0,mask);
    

    
    figure
    
    while length(time_arr)<2^(m-1)*(T_multiples*N_steps)+1
        u = uold - dt/h*.5*(uold.^2 - circshift(uold,1).^2);
        
        
        L1L1R_arr = zeros(size(x));
        L2L2R_arr = zeros(size(x));
        for iq = 1:nq
            tq(iq) = 0.5*dt*xq(iq) + dt/2 +t;
%             [R1L1_iq, R1L1_iq_arr] = compute_temp_1_space_weno(x,dist_x_pl,uold,u,tq(iq),t,dt,mask);
            [R1L1_iq, R1L1_iq_arr, R2L2_iq, R2L2_iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old, max_IUx] = compute_temp_1_space_1(x,dist_x_pl,uold,u,tq(iq),t,dt,mask);

            L1L1R = L1L1R + dt*wq(iq)*R1L1_iq;
            L1L1R_arr = L1L1R_arr + dt*wq(iq)*R1L1_iq_arr;
            L2L2R = L2L2R + dt*wq(iq)*R2L2_iq;
            L2L2R_arr = L2L2R_arr + dt*wq(iq)*R2L2_iq_arr;
            exponential_factor = exponential_factor + wq(iq)*dt*(max_IUx +1);
        end
        
        if it ==0
            c_0_coeff_arr_0 = c_0_coeff_arr_old;
            [error_0, error_0_l2] = space_int_vector_gs_l1(x,dist_x_pl, dist_x_min,u, 0,ex,c_0_coeff_arr_old,mask);
            error_0_l2 = sqrt(error_0_l2);
        end
        
        etat = etat + sum(abs(u - uold))*dt*h;
        etac = etac + 2*(h+dt)*dt*sum(abs(0.5*circshift(uold,1).^2 - 0.5*circshift(uold,0).^2))...
            + (h+dt)^2*dt*h*(sum(abs(circshift(uold,0))+ abs(circshift(uold,1))));
   
        if (showplot==1 && mod(it,fps)==0)
            
            plot(x,u,'b',x,ex(x,t+dt),'r')
            ylim([-.1 .1])
            xlim([-2 2])
            pause(0.01)
        end

        % L1 error with weno reconstruction
%         error_current = space_int_vector_rec_weno_l1(x,dist_x_pl, dist_x_min,u, t,ex);


        % L1 
        if it ==0
            [error_current, error_current_l2] = space_int_vector_gs_l1(x,dist_x_pl, dist_x_min,u, t+dt,ex,c_0_coeff_arr_new,mask);
            
            error_0_l2 = sqrt(error_0_l2);
            bound_arr(1) = error_0;
            bound_ohl_arr(1) = error_0;
            bound_arr_l2(1) = (error_0_l2);
            
%             time_arr(1) = 0;
            error_arr(1) = error_0;
            error_arr_l1l1(1) = error_0;
            error_arr_l2(1) = (error_0_l2);
            EI_arr(1) = 1;
            EI_arr_ohl(1) = 1;
            max_error = error_0;
            max_error_l2 = (error_0_l2);
            
            
            [error_current, error_current_l2] = space_int_vector_gs_l1(x,dist_x_pl, dist_x_min,u, dt,ex,c_0_coeff_arr_new,mask);
            
            error_current_l1l1 = space_int_vector_gs_l1l1(x,dist_x_pl, dist_x_min,uold,u, t,ex,0,mask);
            max_error = max(max_error, error_current);
            max_error_l2 = max(max_error_l2, sqrt(error_current_l2));
            bound_arr(end+1) = L1L1R + error_0;
            bound_arr_l2(end+1) = sqrt(exp(exponential_factor)*(L2L2R+error_0_l2^2));
            bound_ohl_arr(end+1) = (t+1)*(error_0) + etat + etac + ...
                2*t * sqrt(etat + 2*etac); %bound
            error_arr(end+1) = max_error;
            error_arr_l2(end+1) = max_error_l2;
            error_arr_l1l1(end+1) = error_arr_l1l1(end)+error_current_l1l1;
            EI_arr(end+1) = bound_arr(end)/error_arr(end);
            EI_arr_ohl(end+1) = bound_ohl_arr(end)/error_arr(end);
            
        else
            [error_current, error_current_l2] = space_int_vector_gs_l1(x,dist_x_pl, dist_x_min,u, t+dt,ex,c_0_coeff_arr_new,mask);
            
            error_current_l1l1 = space_int_vector_gs_l1l1(x,dist_x_pl, dist_x_min,uold,u, t,ex,0,mask);
            max_error = max(max_error, error_current);
            max_error_l2 = max(max_error_l2, sqrt(error_current_l2));
            bound_arr(end+1) = L1L1R + error_0;
            bound_arr_l2(end+1) = sqrt(exp(exponential_factor)*(L2L2R+error_0_l2^2));
            bound_ohl_arr(end+1) = (t+1)*(error_0) + etat + etac + ...
                2*t * sqrt(etat + 2*etac); %bound
            error_arr(end+1) = max_error;
            error_arr_l2(end+1) = max_error_l2;
            error_arr_l1l1(end+1) = error_arr_l1l1(end)+error_current_l1l1;
            EI_arr(end+1) = bound_arr(end)/error_arr(end);
            EI_arr_ohl(end+1) = bound_ohl_arr(end)/error_arr(end);
        end
        
        t = t+dt;
        it = it+1;
        uold = u;

        time_arr(end+1) = t;
    end
    
    cell_arr{m} = [time_arr;bound_arr;error_arr;EI_arr;bound_ohl_arr;EI_arr_ohl;error_arr_l1l1; error_arr_l2; bound_arr_l2];
end
    save(['FTBS_cell_arr_file_',init_cond,'_burger_try_interp.mat'],'cell_arr')
    save(['FTBS_cell_arr_file_',init_cond,'_burger_try_interp_ul_',num2str(ul),'_ur_',num2str(ur),'_ohlberger_bound.mat'],'cell_arr')

    t_0=0;
    
    fn_plot_FTBS_burger_L1L1_bound(scheme_str,init_cond,t_0,ul,ur)
    fn_plot_FTBS_burger_L1L1_ohlberger_bound(scheme_str,init_cond,t_0,ul,ur)
    fn_plot_FTBS_burger_L1L1_ohlberger_bound_L1time(scheme_str,init_cond,t_0,ul,ur)



function [] = gauss
% For n point gauss quadrature, return evaluation points and weights for
% gauss quadrature over [-1,1]
global nq xq wq
if nq == 1
    xq = 0;
    wq = 1;
elseif nq == 2
    xq = [-0.57735026918962576451, 0.57735026918962576451];
    wq = [0.5, 0.5];
elseif nq == 3
    xq = [-0.7745966692414834, 0, 0.7745966692414834];
    wq = 0.5*[0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
elseif nq > 3
    fprintf(1,'No Gauss quadrature implemented of this degree\n')
    return
end
end