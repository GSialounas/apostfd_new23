function var_out = fn_FTBS_burger(scheme_str,init_conds,ex,t_0)
global nq xq wq % gauss quadrature

showplot = 1;




% N_pts_burger = 100;
% ex = @ (x,t) fn_burger_ex_lin(x,t);%burger_sol(x,t,N_pts_burger);
maxit_arr = [1:3]; % max refinement iterations
% t_0 = 1e-100;


fps = 200;

U = 1;
K = 0;
% T = .2;



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


scheme_arr= {scheme_str};
for i_scheme = 1:length(scheme_arr)
    
    cell_cell_arr_burger = {};
    cell_cell_arr_burger_l1 = {};

    
    
    
    scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
    IC_string = 'sin';
    
    for m = 1:length(maxit_arr)
        Lx= 1; ratio_cTf = 1; part_fine = 1;
        h = 2^(-(maxit_arr(m)+9)); % mesh size normal is +8
        x = create_grid(Lx, ratio_cTf, h, part_fine);
        x2 = -pi +2*pi*(x-x(1))/(x(end)-x(1));
        x=x2(1:end-1);
        dx_fine = x(2)-x(1);
        
        dt = .1*dx_fine^1;
        if m ==1
            if init_conds == 'hatIC'
                T= t_0 + 33*80*dt; % normal is 11*80
                fprintf("The final time T=%5.2f\n",T);
            elseif init_conds == 'sinIC'
                T= t_0 + 33*80*dt; % normal is 11*80
                fprintf("The final time T=%5.2f\n",T);
            elseif init_conds == 'cmbIC'
                T= t_0 + 10*80*dt; % normal is 11*80
                fprintf("The final time T=%5.2f\n",T);
            else
            end
            
        end
        %                 x=x(1:end-1);
        Lx= pi;
        dist_x_pl = circshift(x,-1)-x;
        dist_x_pl(end)= Lx-x(end);
        
        dist_x_min = x- circshift(x,1);
        dist_x_min(1)= Lx-x(end);
        
%         uold = -sin(x);%ex(x,0); % set initial condition
%         u= -sin(x);%ex(x,0); % set initial condition
        
        uold = ex(x,t_0);
        u = uold;
        t = t_0; %initialise time
        L2L2R = 0; %L2L2 accumulation of ||R||^2
        L2L2R_l1 = 0;
        exponential_factor = 0;
        
        % arrays used for EOC for bound and error
        bound_arr = [];%zeros(1,ceil(T/dt));
        l1_bound_arr = [];
        exp_factor_arr = [];
        residual_arr = [];
        error_arr = [];% zeros(1,ceil(T/dt));
        error_arr_l1 = [];
        time_arr = [t_0];% zeros(1,ceil(T/dt));
        dofs_arr = [];% zeros(1,ceil(T/dt));
        EI_index = [];
        const_arr = [];
        error_arr_l1l1 = [];
        bound_arr_ohl = [];
        
        it = 0;
        dx= x(2)-x(1);
        
        etat = 0;
        etac = 0;
        bound_ohl = 0;
        error_0 = 0;
        
        
        mask = ones(size(x));
        if init_conds == 'stpIC'
            mask(abs(x)>1) =0;
        elseif init_conds == 'hatIC'
        elseif init_conds == 'cmbIC'
            mask(abs(x)>2.0) = 0;
        else
        end
        
        while t<=T-dt/2
            
            %             spat_disc = 'BS';
            %             f_h_old = flux_central(dist_x_pl,dt,uold);%flux_fn(dist_x_pl,dt,uold);
            %             u_stage_1 = uold - dt*uold.*f_h_old;
            %             u = .5 * uold+ .5 * u_stage_1  - .5 * dt * u_stage_1.* flux_central(dist_x_pl,dt,u_stage_1);
            %             f_h_new = flux_central(dist_x_pl,dt,u);%;flux_fn(dist_x_pl,dt,u);%1/dx_fine*(flux_eo(u,circshift(u,-1),dx_fine)- flux_eo(circshift(u,1),circshift(u,0),dx_fine));
            uold= u;
            spat_disc = "LxF";
            c=1;
            if spat_disc =="WENO5"
                f_h_old  = WENO5resAdv1d_fdm_gs(uold,flux,dflux,S,dx);%resWENO5(uold_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu); %WENO5resAdv1d_fdm_gs(uold_nu,flux,dflux,S,dx);%
                
                % Three stage
                ustage_0 = uold;
                
                ustage_1 = uold - dt*f_h_old;
                f_h_stage_1 =  WENO5resAdv1d_fdm_gs(ustage_1,flux,dflux,S,dx);%resWENO5(ustage_1,c,dx);% fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_1); %WENO5resAdv1d_fdm_gs(ustage_1,flux,dflux,S,dx);%
                
                ustage_2 = (.75)*uold +.25*ustage_1 - .25 *dt *(f_h_stage_1);
                f_h_stage_2 = WENO5resAdv1d_fdm_gs(ustage_2,flux,dflux,S,dx);%resWENO5(ustage_2,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_2); %WENO5resAdv1d_fdm_gs(ustage_2,flux,dflux,S,dx);%
                
                ustage_3 = (1/3)*uold +(2/3)*ustage_2 - (2/3)*dt*(f_h_stage_2);
                
                u = ustage_3;
                [f_h_new, f_u] = WENO5resAdv1d_fdm_gs(u,flux,dflux,S,dx);%resWENO5(u_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_3); %WENO5resAdv1d_fdm_gs(u_nu,flux,dflux,S,dx);%
            elseif spat_disc== "WENO3"
                f_h_old  = WENO3resAdv1d(uold,flux,dflux,S,dx);%resWENO5(uold_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu); %WENO5resAdv1d_fdm_gs(uold_nu,flux,dflux,S,dx);%
                
                % Three stage
                ustage_0 = uold;
                
                ustage_1 = uold - dt*f_h_old;
                f_h_stage_1 =  WENO3resAdv1d(ustage_1,flux,dflux,S,dx);%resWENO5(ustage_1,c,dx);% fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_1); %WENO5resAdv1d_fdm_gs(ustage_1,flux,dflux,S,dx);%
                
                ustage_2 = (.75)*uold +.25*ustage_1 - .25 *dt *(f_h_stage_1);
                f_h_stage_2 = WENO3resAdv1d(ustage_2,flux,dflux,S,dx);%resWENO5(ustage_2,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_2); %WENO5resAdv1d_fdm_gs(ustage_2,flux,dflux,S,dx);%
                
                ustage_3 = (1/3)*uold +(2/3)*ustage_2 - (2/3)*dt*(f_h_stage_2);
                
                u = ustage_3;
                [f_h_new,f_u] = WENO3resAdv1d(u,flux,dflux,S,dx);%resWENO5(u_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_3); %WENO5resAdv1d_fdm_gs(u_nu,flux,dflux,S,dx);%
%                 
%                 f_h_old = (.5/dx)*(uold.^2-circshift(uold,1).^2);
%                 u = uold-dt*f_h_old;
%                 f_h_new= (.5/dx)*(u.^2-circshift(u,1).^2);
            elseif spat_disc == "LxF"
                f_h_old = flux_lxf(dist_x_pl,dt,uold);
                u = uold - dt*f_h_old;
                f_h_new =  flux_lxf(dist_x_pl,dt,u);
            elseif spat_disc == "FTBS"
                f_h_old = (1/dx)*.5*(uold.^2-circshift(uold,1).^2);%flux_lxf(dist_x_pl,dt,uold);
                u = uold - dt*f_h_old;
                f_h_new =  (1/dx)*.5*(u.^2-circshift(u,1).^2);%flux_lxf(dist_x_pl,dt,u);
            else
            end
            
            
            % Initial indicator calculation at time step t = it*dt
            L2L2R_arr = zeros(size(x));
            L2L2R_l1_arr = zeros(size(x));
            L2L2R_arr_t = zeros(size(x));
            L2L2R_arr_x = zeros(size(x));
            for iq = 1 : nq
                tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
%                 [RL2iq,RL2iq_arr, RL2iq_l1,RL2iq_l1_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,R_h_burger,max_IUx] = compute_burger_rec_space_3_time_3_alt_rec(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, f_h_old, f_h_new,spat_disc,c); %compute R at gauss point

                [RL2iq,RL2iq_arr, RL2iq_l1,RL2iq_l1_arr, c_0_coeff_arr_new, c_0_coeff_arr_old, max_IUx,xiq,IUr,IUt,IUx] = compute_Rs_vector_temp_1_spatiotemp_1_burger(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, f_h_old, f_h_new,spat_disc,c,ex); %compute R at gauss point
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
            
            etat = etat + sum(abs(u - uold))*dt*h;
            etac = etac + 2*(h+dt)*dt*sum(abs(0.5*circshift(uold,1).^2 - 0.5*circshift(uold,0).^2))...
                + (h+dt)^2*dt*h*(sum(abs(circshift(uold,0))+ abs(circshift(uold,1))));
 
            
            if it ==0
                c_0_coeff_arr_0 = c_0_coeff_arr_old;
                x_0 = x;
                dist_x_pl_0 = dist_x_pl;
                dist_x_min_0 = dist_x_min;
%                 error_0 =space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex);%%( space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0));%
%                 error_0_l1 =space_int_vector_rec_weno_burger_l1(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex);%%( space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0));%
                error_0 = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
                error_0_l1 =(space_int_vector_gs_l1(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
%             error_0_l1 = 0;
            end
            
            if (showplot && mod(it-1,fps)==0)
                l_coef= length(c_0_coeff_arr_new);
%                 IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl;dist_x_pl.^2;dist_x_pl.^3],1);
                IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl],1);
%                 subplot(1,3,1)
                plot(x,u,'b',x,ex(x,t+dt),'r') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                axis([x(1) x(end) -1 1])
%                 subplot(1,3,2)
%                 plot(xiq,IUt,'b',xiq,IUr.*IUx,'r')
%                 axis([x(1) x(end) -2 2])
% %                 plot(x,f_h_new,'r*',x,.5*(ex(x,t).^2),'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
%                 subplot(1,3,3)
%                 plot(xiq,1/dt*(u-uold),'r')
                pause(0.01)
            end
            
            if it ==0
                % t = 0
                const_arr(1) = 1;
                
%                 bound_arr(1) = sqrt(space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex));%sqrt( space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0));%
                bound_arr(1) = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
                exp_factor_arr(1)= sqrt(exp(exponential_factor));
                residual_arr(1) = bound_arr(1);
%                 error_arr(1) = sqrt(space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex));%sqrt( space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0)) ;%
                error_arr(1) =  sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
                error_arr_l1(1) = (space_int_vector_gs_l1(x,dist_x_pl,dist_x_min,uold,t_0,ex,c_0_coeff_arr_0));
%                 error_arr_l1(1) = 0;
                l1_bound_arr(1) = error_0_l1;
                
                error_old = max(error_arr);
                EI_index(1) = bound_arr(1)/error_arr(1);
                bound_arr_ohl(1) = error_0_l1;
                error_arr_l1l1(1) = error_0_l1;
                
                
                % we will calculate the quantities for the L1 bound here
                % We need four quantities: eta_0, eta_t, eta_c and
                % eta_overbar
                
                % eta_0 is essentially 
%                 eta_0 = space_int_vector_rec_weno_l1(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex);
%                 eta_overbar = 0; % this is the boundary condition
%                 eta_t = eta_t + h*dt*sum(abs(u-uold));
%                 eta_c = eta_c + 2*(h+dt)*dt*(max(f_h_new)); 
%                 eta_c = eta_c + 2*(h+dt)*dt*(1/h)*(max(.5*(u.^2-circshift(u.^2,1))))+ 2*(h+dt)*dt*(1/h)*(max(.5*(-u.^2+circshift(u.^2,1)))); 
%                 eta_c = eta_c + 2*(h+dt)*dt*(max(.5*(-u.^2+circshift(u.^2,1)))); 
                
                %2*T*1*(eta_0 );
                
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
                bound_arr_ohl(end+1) = error_0_l1 + (T)*(error_0_l1) + etat + etac + ...
                2*T * sqrt(1/(T) +2)* sqrt(etat + 2*etac); %bound
                error_current_l1l1 = space_int_vector_gs_l1l1(x,dist_x_pl, dist_x_min,uold,u, t,ex,0,mask);
                error_arr_l1l1(end+1) = error_arr_l1l1(end) + error_current_l1l1;
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
                
%                 eta_0 = space_int_vector_rec_weno_l1(x,dist_x_pl,dist_x_min,uold,(1e-100)*dt,ex);
%                 eta_overbar = 0; % this is the boundary condition
%                 eta_t = eta_t + h*dt*sum(abs(u-uold));
%                 eta_c = eta_c + 2*(h+dt)*dt*(max(f_h_new)); 
%                 eta_c = eta_c + 2*(h+dt)*dt*(1/h)*(max(.5*(u.^2-circshift(u.^2,1)))) + 2*(h+dt)*dt*(1/h)*(max(.5*(-u.^2+circshift(u.^2,1))));
%                 eta_c = eta_c + 2*(h+dt)*dt*(max(.5*(-u.^2+circshift(u.^2,1))));
                
                l1_bound_arr(end+1) =error_0_l1 + L2L2R_l1;%2*T*1*(eta_0 + eta_overbar + (eta_t+eta_c)*1 + eta_c*1 );
                bound_arr_ohl(end+1) = error_0_l1 + ((T))*(error_0_l1) + etat + etac + ...
                    2*(T) * sqrt(1/(T) +2) * sqrt(etat + 2*etac); %bound
                error_current_l1l1 = space_int_vector_gs_l1l1(x,dist_x_pl, dist_x_min,uold,u, t,ex,0,mask);
                error_arr_l1l1(end+1) = error_arr_l1l1(end) + error_current_l1l1;
                

            end
            
            t = t +dt;
            it=it+1;
%             uold = u;
%             max(L2L2R_arr_x);
            time_arr(end+1) = t_0+ it*dt;
            
%             t = t +dt;
%             uold = u;
%             max(L2L2R_arr_x);
%             time_arr(it) = it*dt;
%             IUx = sum(c_0_coeff_arr_new(2:end, :).*[ones(size(x));2*dist_x_pl; 3*dist_x_pl.^2],1);
%             
% %             bound_arr(it) = sqrt(exp(it*dt*(1+max(abs(IUx))))*(L2L2R + space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0)));
%             bound_arr(it) = sqrt(exp(exponential_factor)*(L2L2R + space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0)));
% 
%             % error_arr_eoc(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,it*dt,ex));
%             error_arr(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,it*dt,ex,c_0_coeff_arr_new));%sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,it*dt,ex,c_0_coeff_arr_new));
%             dofs_arr(it) = length(x);
%             EI_index(it) = bound_arr(it)./error_arr(it);
%             it = it + 1;
        end
%         finalL2err(m) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,(it-1)*dt,ex,c_0_coeff_arr_new)); %compute final time L2 error
% %         R(m) = sqrt(exp((it-1)*dt*(1+max(abs(IUx))))*(space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex,10^-10,ex,c_0_coeff_arr_0) + L2L2R)); %bound
%         R(m) = sqrt(exp(exponential_factor)*(space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex,10^-10,ex,c_0_coeff_arr_0) + L2L2R)); %bound
%         finalL2err(m) = space_int_vector_rec_weno_burger_l1(x,dist_x_pl,dist_x_min,u,(it)*dt,ex);%sqrt( space_int_vector_rec_weno(x,dist_x_pl,dist_x_min,u,(it)*dt,ex));% sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,(it)*dt,ex,c_0_coeff_arr_new)); %
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
        fprintf(1,'      ||R||_L2(L1) = %.5f  EOC = %1.2f\n',R(m),EOCR(m))
        fprintf(1,'                EI = %.5f  \n',EI(m))
        
        % Arrays for plotting
        cell_cell_arr_burger{m}=[time_arr;bound_arr;error_arr;EI_index; exp_factor_arr; residual_arr;const_arr];
        cell_cell_arr_burger_l1{m}=[time_arr;bound_arr;error_arr;l1_bound_arr;error_arr_l1;EI_index; exp_factor_arr; residual_arr;const_arr];
        cell_cell_arr_burger_l1l1{m}=[time_arr;bound_arr;error_arr;l1_bound_arr;error_arr_l1;bound_arr_ohl;error_arr_l1l1 ];

    end
    
    
    save([scheme_str,'_cell_arr_file_sinIC_burger.mat'],'cell_cell_arr_burger')
    save([scheme_str,'_cell_arr_file_sinIC_burger_l1_bound.mat'],'cell_cell_arr_burger_l1')
    save([scheme_str,'_cell_arr_file_',init_conds,'_burger_l1l2_bound_comparison.mat'],'cell_cell_arr_burger_l1')
    save([scheme_str,'_cell_arr_file_',init_conds,'_burger_l1l1_bound_comparison.mat'],'cell_cell_arr_burger_l1l1')

    
end
    var_out=1;
end



