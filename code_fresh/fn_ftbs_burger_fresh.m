function var_out = fn_ftbs_burger_fresh(scheme_str,init_conds,ex,t_0)
global nq xq wq

showplot=1;
fps=100;

maxit_arr =  1:3;
L=1;
c=1;
spat_disc= "BS";
fluxfun='burgers'; % select flux function
% Define our Flux function
switch fluxfun
    case 'linear'   % Scalar Advection, CFL_max: 0.65
        c=1; flux = @(w) c*w;
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.40
        flux = @(w) w.^2/2;
        dflux = @(w) w;
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

cell_cell_arr_burger = {};
cell_cell_arr_burger_l1 = {};
   

for m  = 1:length(maxit_arr)
    h = 2^(-(maxit_arr(m)+8));
    x = -1:h:1;
    x2 = -pi+2*pi*(x-x(1)/(x(end)-x(1)));
    x = x2;%(1:end-1);
    dist_x_pl = x(2:end)-x(1:end-1);
    dist_x_min = dist_x_pl;
%     x=x(1:end-1);
    dx = x(2)-x(1);
    dt = .1*(x(2)-x(1));
    if m==1
        T = 11*80*dt;
    end
    
    
    Lx = pi;
        x= x(1:end-1);

    uold = ex(x,0);
    u = uold;
    t = 0;
    L2L2R = 0; %L2L2 accumulation of ||R||^2
    L2L2R_l1 = 0;
    exponential_factor = 0;
    bound_arr = [];
    l1_bound_arr =[];
    
    exp_factor_arr = [];
    residual_arr = [];
    error_arr = [];% zeros(1,ceil(T/dt));
    error_arr_l1 = [];
    time_arr = [t_0];% zeros(1,ceil(T/dt));
    dofs_arr = [];% zeros(1,ceil(T/dt));
    EI_index = [];
    const_arr = [];
    it = 0;
    
    while t<= T-dt/2
        uold = u;
        f_h_old = (1/dx)*.5*(uold.^2-circshift(uold,1).^2);
        f_h_old = (1/dx)*.5*(uold(2:end).^2-uold(1:end-1).^2);

        u(2:end) = uold(2:end) - dt*f_h_old;
        f_h_new = (1/dx)*.5*(u(2:end).^2-u(1:end-1).^2);
        
%         f_h_old = (1/dx)*1*(uold(2:end)-uold(1:end-1));
% 
%         u(2:end) = uold(2:end) - dt*f_h_old;
%         f_h_new = (1/dx)*1*(u(2:end)-u(1:end-1));
        
        % Initial indicator calculation at time step t = it*dt
        L2L2R_arr = zeros(size(x));
        L2L2R_l1_arr = zeros(size(x));
        L2L2R_arr_t = zeros(size(x));
        L2L2R_arr_x = zeros(size(x));
        
        for iq = 1 : nq
            tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
            %                 [RL2iq,RL2iq_arr, RL2iq_l1,RL2iq_l1_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,R_h_burger,max_IUx] = compute_burger_rec_space_3_time_3_alt_rec(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, f_h_old, f_h_new,spat_disc,c); %compute R at gauss point
            
            
            [RL2iq,RL2iq_arr, RL2iq_l1,RL2iq_l1_arr, c_0_coeff_arr_new, c_0_coeff_arr_old, max_IUx,xiq,IUr,IUt,IUx] = compute_Rs_vector_temp_1_spatiotemp_1_burger(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, f_h_old, f_h_new,spat_disc,flux,ex); %compute R at gauss point

%             [RL2iq,RL2iq_arr, RL2iq_l1,RL2iq_l1_arr, c_0_coeff_arr_new, c_0_coeff_arr_old, max_IUx,xiq,IUr,IUt,IUx] = compute_Rs_vector_temp_1_spatiotemp_1_burger_fresh_non_periodic(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, f_h_old, f_h_new,spat_disc,c,ex); %compute R at gauss point
            L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
            L2L2R_arr = L2L2R_arr +wq(iq)*dt*(RL2iq_arr);
            
            L2L2R_l1 = L2L2R_l1 + wq(iq)*dt*(RL2iq_l1); %quadrature formula
            L2L2R_l1_arr = L2L2R_l1_arr +wq(iq)*dt*(RL2iq_l1_arr);
            % We will now also make a calculation of the factor that
            % goes into the exponential
            exponential_factor = exponential_factor + wq(iq)*dt*(max_IUx +1); %quadrature formula
            
            
            if (L2L2R<0)
                disp('L2L2R<0')
            end
        end
        
        if it ==0
            c_0_coeff_arr_0 = c_0_coeff_arr_old;
            x_0 = x;
            dist_x_pl_0 = dist_x_pl;
            dist_x_min_0 = dist_x_min;
            %                 error_0 =space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex);%%( space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0));%
            %                 error_0_l1 =space_int_vector_rec_weno_burger_l1(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex);%%( space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0));%
            
            % the error_0_l2 is not square rooted
%             [error_0_l1, error_0_l2] = space_int_vector_gs_l1_non_periodic(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0);
%             error_0_l2 = sqrt(error_0_l2);
            error_0 = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
            error_0_l1 =(space_int_vector_gs_l1(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
        end
        
        
        if (showplot && mod(it-1,fps)==0)
            l_coef= length(c_0_coeff_arr_new);
            %                 IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl;dist_x_pl.^2;dist_x_pl.^3],1);
            IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl],1);
%             subplot(1,3,1)
            plot(x,IU,'r',x,ex(x,t+dt),'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
%             axis([0 0.1 -.1 1.1])
            legend('IU','exact')
            title('exact vs IU')
            xlabel('x')
            ylabel('solution')
%             subplot(1,3,2)
%             plot(xiq,IUt,'b',xiq,IUr.*IUx,'r')
%             axis([x(1) x(end) -2 2])
%             %                 plot(x,f_h_new,'r*',x,.5*(ex(x,t).^2),'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
%             subplot(1,3,3)
%             plot(xiq,1/dt*(u(1:end-1)-uold(1:end-1)),'r')
            pause(0.01)
        end
        
        if it ==0
            % t = 0
            const_arr(1) = 1;
            
            %                 bound_arr(1) = sqrt(space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex));%sqrt( space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0));%
            bound_arr(1) =error_0;%sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
            exp_factor_arr(1)= sqrt(exp(exponential_factor));
            residual_arr(1) = bound_arr(1);
            %                 error_arr(1) = sqrt(space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex));%sqrt( space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0)) ;%
            error_arr(1) =  error_0;%sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
            error_arr_l1(1) = error_0_l1;%(space_int_vector_gs_l1(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
            
            error_old = max(error_arr);
            EI_index(1) = bound_arr(1)/error_arr(1);
            
            l1_bound_arr(1) = error_0_l1;%2*T*1*(eta_0 );
            
            % t = dt
            const_arr(end+1) =1;
            bound_arr(end+1)= sqrt(exp(exponential_factor)*(L2L2R + error_0^2));%sqrt(exp(exponential_factor)*(L2L2R + space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0)));%
            exp_factor_arr(end+1)= sqrt(exp(exponential_factor));
            residual_arr(end+1) = sqrt((L2L2R + error_0));%sqrt(L2L2R + space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0));
            %                 error_arr(end+1) = sqrt( space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,u,1*dt,ex));% sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,1*dt,ex,c_0_coeff_arr_new));%
            error_arr(end+1) = max(error_arr(1),sqrt( space_int_vector_gs(x,dist_x_pl,dist_x_min,u,t_0+ dt,ex,c_0_coeff_arr_new)));
            error_arr_l1(end+1) = max(error_arr_l1(1),(space_int_vector_gs_l1(x,dist_x_pl,dist_x_min,u,t_0+ dt,ex,c_0_coeff_arr_new)));
            %                 error_arr_l1(end+1) = space_int_vector_rec_weno_burger_l1(x,dist_x_pl,dist_x_min,u,1*dt,ex);
            error_old = max(error_arr);
            EI_index(end+1) = bound_arr(end)./error_arr(end);
            l1_bound_arr(end+1) = error_0_l1 + L2L2R_l1;%2*T*1*(eta_0 + eta_overbar + (eta_t+eta_c)*1 + eta_c*1 );
        else
            const_arr(end+1) = 1;
            bound_arr(end+1) = sqrt(exp(exponential_factor)*(L2L2R + error_0^2));% sqrt(exp(exponential_factor)*(L2L2R + space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0)));%
            exp_factor_arr(end+1)= sqrt(exp(exponential_factor));
            residual_arr(end+1) = sqrt(L2L2R + error_0);
            % error_arr(end+1) = sqrt( space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,u,(it+1)*dt,ex)); % sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,(it+1)*dt,ex,c_0_coeff_arr_new));%
            % error_arr_l1(end+1) = sqrt( space_int_vector_rec_weno_burger_l1(x,dist_x_pl,dist_x_min,u,(it+1)*dt,ex));
            error_max = max(error_arr);
            error_arr(end+1) = max(error_max,sqrt( space_int_vector_gs(x,dist_x_pl,dist_x_min,u, t_0+(it+1)*dt,ex,c_0_coeff_arr_new)));
            error_l1_max = max(error_arr_l1);
            error_arr_l1(end+1) = max(error_l1_max,(space_int_vector_gs_l1(x,dist_x_pl,dist_x_min,u,t_0+(it+1)*dt,ex,c_0_coeff_arr_new)));
            
            error_old = max(error_arr);
            EI_index(end+1) = bound_arr(end)./error_arr(end);
            l1_bound_arr(end+1) =error_0_l1 + L2L2R_l1;%2*T*1*(eta_0 + eta_overbar + (eta_t+eta_c)*1 + eta_c*1 );
            
        end
        
        t=t+dt;
        it = it+1;
        time_arr(end+1) = t_0+it*dt;
    end
    finalL2err(m) = space_int_vector_gs_l1(x,dist_x_pl,dist_x_min,u,t_0+(it)*dt,ex,c_0_coeff_arr_new);
    %         R(m) = sqrt(exp((it-1)*dt*(1+max(abs(IUx))))*(space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex,10^-10,ex,c_0_coeff_arr_0) + L2L2R)); %bound
    R(m) = error_0_l1 + L2L2R_l1;%sqrt(exp(exponential_factor)*(error_0 + L2L2R));%sqrt(exp(exponential_factor)*(space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex,10^-10,ex,c_0_coeff_arr_0) + L2L2R)); %
    
    EOCe(1) = 0;
    EOCR(1) = 0;
    if m > 1
        EOCe(m) = log(finalL2err(m-1)/finalL2err(m))/log(2);
        EOCR(m) = log(R(m-1)/R(m))/log(2);
        EOC_error_arr(m) = EOCe(m);
        EOC_bound_arr(m) = EOCR(m);
    end
    EI(m) = R(m)/finalL2err(m);
    fprintf(1,'||(u - IU)(T)||_L1 = %.5f  EOC = %1.2f\n',finalL2err(m),EOCe(m))
    fprintf(1,'      ||R||_L1(L1) = %.5f  EOC = %1.2f\n',R(m),EOCR(m))
    fprintf(1,'                EI = %.5f  \n',EI(m))
    
    % Arrays for plotting
    cell_cell_arr_burger{m}=[time_arr;bound_arr;error_arr;EI_index; exp_factor_arr; residual_arr;const_arr];
    cell_cell_arr_burger_l1{m}=[time_arr;bound_arr;error_arr;l1_bound_arr;error_arr_l1;EI_index; exp_factor_arr; residual_arr;const_arr];
    
    
    
end
 
    save([scheme_str,'_cell_arr_file_sinIC_burger_fresh.mat'],'cell_cell_arr_burger')
    save([scheme_str,'_cell_arr_file_sinIC_burger_l1_bound_fresh.mat'],'cell_cell_arr_burger_l1')
    save([scheme_str,'_cell_arr_file_',init_conds,'_burger_l1l2_bound_comparison_fresh.mat'],'cell_cell_arr_burger_l1')

var_out=1;
end