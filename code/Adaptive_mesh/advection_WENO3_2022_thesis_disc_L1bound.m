% This script is dated 30th May 2022 and it is a test that will go in the
% phd thesis.  We will use this to test the post-shock regime L1 bound
% using advection with discontinuous IC.


clear all;
close all;
clc;
mesh_adapt=0;
if mesh_adapt == 1
    adaptivity = 'adapt_ON';
else
    adaptivity = 'adapt_OFF';
end
maxit_arr = [1:4]; % max refinement iterations
record_video_arr = zeros(size(maxit_arr));
showplot_arr = zeros(size(maxit_arr));

record_video_arr(end)= 0;
showplot_arr(end) =0;
% record_video = 0;
fps= 20;

rep_video = 1;
f_size = 16;
max_number_ref_levels = 1;
gamma_param = .5;

scheme_arr = {'WENO3_nonu'};
intrpl_arr = {'P3'};
init_conds = 'smooth';
% interpolant_arr
exponent_arr = [1];
global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();
showplot=1;



fluxfun='linear'; % select flux function
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

for i_scheme = 1:length(scheme_arr)
    
    
    cell_cell_arr = {};
    figgy_parasite = figure();
    
    
    scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
    interpolant = intrpl_arr{1};
    exponent =exponent_arr(1);
    IC_string = 'sin';
    
    T = .5; % final time
    
    Lx = 1;
    ratio_cTf = 1;
    part_fine = 1;
    
    %     ex = @(x,t) sin(2*pi*(x-t)); %smooth exact solution
    %     ex_x =  @(x,t) 2*pi*cos(2*pi*(x-t));
    
    
    cntr = .25;
    Lx=1;
    pow_exp=100;
    %             ratio_cTf=2;
    part_fine=1;
    radius = 0.125;
    % ex = @(x,t) exp(-100*(x-cntr - t).^2); %smooth exact solution
    
    if init_conds == 'stepIC'
        ex = @(x,t) fn_stp_exact(x,cntr,t,Lx,radius) ;%
    elseif init_conds == 'smooth'
        ex = @(x,t) sin(2*pi*(x-t)); %smooth exact solution
    elseif init_conds =='hatfIC'
        h_kmin = .25;
        h_kpl = .25;
        ex = @(x,t) fn_hat_exact(x,cntr,t,Lx,h_kmin,h_kpl);
    else
    end
    %  ex = @(x,t) exp(-10*((10*(x-cntr - t))).^2); %smooth exact solution
    %  grid_x_sin=@(x) (sin(50*pi*x)).^2;
    %   ex = @(x,t) max(exp(-pow_exp*((mod(x-cntr-t,Lx)).^2)),exp(-pow_exp*((mod(x-cntr-t,-Lx)).^2))); %smooth exact solution
    
    
    fprintf(1,'Using the %s scheme\n',scheme)
    maxit = length(maxit_arr);
    
    
    for m =1:length(maxit_arr)
        record_video = 0;
        %                 showplot =showplot_arr(m);
        
        
        % produce initial mesh
        h = 2^(-(maxit_arr(m)+4)); % mesh size
        %         dt_coarsest = .1 * 2^(-(maxit_arr(end)+8)); % mesh size
%         if mesh_adapt ==0
%             x= linspace(0,1,741);
%             
%         else
            x = create_grid(Lx, ratio_cTf, h, part_fine);
            
%         end
        dt =  .1 * h;%.1*h^1; % initial time-step
        dx_arr(m) = h;
        dt_arr(m) = dt; % timestep size
        
        dist_x_pl = x(2:end)-x(1:end-1);
        dist_x_min = circshift(dist_x_pl,1);
        x=x(1:end-1);

 
        
        
        
        uold = ex(x,0); % set initial condition
        
        t = 0; %initialise time
        it = 0;
        i = 0;
        L2L2R = 0; %L2L2 accumulation of ||R||^2
        L2L2R_arr = zeros(size(x)); %element-wise L2L2 accumulation of ||R||^2: i.e. ||R||_(element)^2
        reps = 0;
        
        % arrays used for EOC for bound and error
        bound_arr = []; %zeros(1,ceil(T/dt));
        error_arr = []; % zeros(1,ceil(T/dt));
        time_arr  = []; % zeros(1,ceil(T/dt));
        dofs_arr  = []; % zeros(1,ceil(T/dt));
        EI_index  = [];
        
        time_arr(1) = 0;
        
        u= uold;
        
        
        %         L2L2_arr_cumulative= zeros(size(x));
        error_old = 0;
        error_new = 0;
        if m ==1
            dt_coarsest=dt;
        end
        error_new_arr=[0];
        
%         % form initial mesh
%         tree_list = [];
%         root_ids = [];
%         x = [x, Lx];
%         tol = 1e-6;
%         max_ref= 4;
%         for i = 1:length(x)-1
%             if i==1
%                 tree_list = [ Tree(element(x(i), x(i+1),i,tol),i,tol)];
%                 %                 tree_list(i) = tree_list(i).assign_vals([sin(x(i)), sin(x(i+1))], exp(-100*(x(i)-.5).^2));
%                 root_ids(i) = i;
%             else
%                 tree_list(end+1) = Tree(element(x(i), x(i+1),i,tol),i,tol);
%                 %                 tree_list(i) = tree_list(i).assign_vals([sin(x(i)), sin(x(i+1))], exp(-100*(x(i)-.5).^2));
%                 root_ids(i) = i;
%             end
%         end
%         
%         if mesh_adapt==1
%             M = mesh(tree_list,root_ids,max_ref);
%             M.updateMesh();
%             x = x(1:end-1);
%             M.asgnVals(uold, zeros(size(uold)));
%             
%             M.refineGlobal();
%             M.refineGlobal();
%             M.refineGlobal();
%             M.refineGlobal();
%             [x_ref,uold_ref] = M.getNewVals();
%             
%             % amend the grid and uld
%             x_ref = [x_ref, Lx];
%             uold  = uold_ref;
%             
%             dist_x_pl_ref = x_ref(2:end) -x_ref(1:end-1);
%             dist_x_min_ref = circshift(dist_x_pl_ref, 1);
%             dist_x_pl = dist_x_pl_ref;
%             dist_x_min = dist_x_min_ref;
%             x = x_ref(1:end-1);
%         else
%             x=x(1:end-1);
%         end
%         if mesh_adapt == 1
%             uold= ex(x,0);
%             M.asgnVals(uold, zeros(size(uold)));
%         end
        
        while t <= T-dt/2          
            % initial solve at time step  t = it*dt
            spatial_disc = "WENO3_nonu";
            
            % Initial solve and initial estimator computation
            [~, uold_x] = WENO3resAdv1d_periodic_nonu(x, x,dist_x_pl, uold,flux,dflux,Lx);
            % Three stage
            ustage_0 = uold;
            ustage_1 = uold - dt*uold_x;
            
            [~, uold_x_stage_1] = WENO3resAdv1d_periodic_nonu(x, x,dist_x_pl, ustage_1,flux,dflux,Lx);
            ustage_2 = (.75)*uold +.25*ustage_1 - .25 *dt *(uold_x_stage_1);
            
            [~, uold_x_stage_2] = WENO3resAdv1d_periodic_nonu(x, x,dist_x_pl, ustage_2,flux,dflux,Lx);
            ustage_3 = (1/3)*uold +(2/3)*ustage_2 - (2/3)*dt*(uold_x_stage_2);
           
            u = ustage_3;
            [~, u_x] = WENO3resAdv1d_periodic_nonu(x, x,dist_x_pl, ustage_3,flux,dflux,Lx);
            % Initial indicator calculation at time step t = it*dt
            
            if mesh_adapt==0
                L2L2R_arr= zeros(size(x));
                L2L2R_arr_t= zeros(size(x));
                L2L2R_arr_x= zeros(size(x));
                for iq = 1 : nq
                    tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                    if interpolant=='P1'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_1(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    elseif interpolant=='P2'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_2(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    elseif interpolant=='P3'
                        [RL2iq, RL2iq_arr,c_0_coeff_arr_new, c_0_coeff_arr_old, IUt, IUx, IUt_exact, IUx_exact] = compute_rec_space_3_time_3_alt_rec_nonu(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc,c,Lx);%compute R at gauss point
                    else
                    end
                    L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                    L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
%                     L2L2R_arr_t = RL2iq_arr_t;%L2L2R_arr_t + wq(iq)*dt*(RL2iq_arr_t); %quadrature formula
%                     L2L2R_arr_x = RL2iq_arr_x;% L2L2R_arr_x + wq(iq)*dt*(RL2iq_arr_x); %quadrature formula
                    
                    if (L2L2R<0)
                        disp('L2L2R<0')
                    end
                end
                
            %% refine the mesh
            else
                L2L2R_arr= zeros(size(x));
                L2L2R_arr_t= zeros(size(x));
                L2L2R_arr_x= zeros(size(x));
                for iq = 1 : nq
                    tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                    if interpolant=='P1'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_1(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    elseif interpolant=='P2'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_2(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    elseif interpolant=='P3'
                        [RL2iq, RL2iq_arr,c_0_coeff_arr_new, c_0_coeff_arr_old, IUt, IUx, IUt_exact, IUx_exact] = compute_rec_space_3_time_3_alt_rec_nonu(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc,c,Lx); %compute R at gauss point

                    else
                    end
%                     L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                    L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
%                     L2L2R_arr_t = RL2iq_arr_t;%L2L2R_arr_t + wq(iq)*dt*(RL2iq_arr_t); %quadrature formula
%                     L2L2R_arr_x = RL2iq_arr_x;% L2L2R_arr_x + wq(iq)*dt*(RL2iq_arr_x); %quadrature formula
                    
                    if (L2L2R<0)
                        disp('L2L2R<0')
                    end
                end
                if it ==0 && mesh_adapt==1
                    c_0_coeff_arr_0 = c_0_coeff_arr_old;
                    x_0 = x;
                    dist_x_pl_0 = dist_x_pl;
                    dist_x_min_0 = dist_x_min;
                    error_0 = space_int_vector_rec_weno_nonu(x,dist_x_pl,dist_x_min,uold,0*dt,ex,Lx);

                end
                
                M.asgnVals(uold, RL2iq_arr);
                M.mark_ref_coar();
                M.ref_coarsen_local();
                [x_ref, uold_ref] = M.getNewVals();
                
                % amend the grid and uld
                x_ref = [x_ref, Lx];
                uold  = uold_ref;
                
                dist_x_pl_ref = x_ref(2:end) -x_ref(1:end-1);
                dist_x_min_ref = circshift(dist_x_pl_ref, 1);
                dist_x_pl = dist_x_pl_ref;
                dist_x_min = dist_x_min_ref;
                x = x_ref(1:end-1);
                
                % take the time-step with the refined uold
                
                %                 dt = .1 *min(dist_x_pl);
                [~, uold_x] = WENO3resAdv1d_periodic_nonu(x, x,dist_x_pl, uold,flux,dflux,Lx);
                % Three stage
                ustage_0 = uold;
                ustage_1 = uold - dt*uold_x;
                
                [~, uold_x_stage_1] = WENO3resAdv1d_periodic_nonu(x, x,dist_x_pl, ustage_1,flux,dflux,Lx);
                ustage_2 = (.75)*uold +.25*ustage_1 - .25 *dt *(uold_x_stage_1);
                
                [~, uold_x_stage_2] = WENO3resAdv1d_periodic_nonu(x, x,dist_x_pl, ustage_2,flux,dflux,Lx);
                ustage_3 = (1/3)*uold +(2/3)*ustage_2 - (2/3)*dt*(uold_x_stage_2);
                
                u = ustage_3;
                [~, u_x] = WENO3resAdv1d_periodic_nonu(x, x,dist_x_pl, ustage_3,flux,dflux,Lx);
                
                % bound for plotting purposes
                
                L2L2R_arr= zeros(size(x));
                L2L2R_arr_t= zeros(size(x));
                L2L2R_arr_x= zeros(size(x));
                for iq = 1 : nq
                    tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                    if interpolant=='P1'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_1(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    elseif interpolant=='P2'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_2(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    elseif interpolant=='P3'
                        [RL2iq, RL2iq_arr,c_0_coeff_arr_new, c_0_coeff_arr_old, IUt, IUx, IUt_exact, IUx_exact] = compute_rec_space_3_time_3_alt_rec_nonu(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc,c,Lx);%compute R at gauss point
                    else
                    end
                    L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                    L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
%                     L2L2R_arr_t = RL2iq_arr_t;%L2L2R_arr_t + wq(iq)*dt*(RL2iq_arr_t); %quadrature formula
%                     L2L2R_arr_x = RL2iq_arr_x;% L2L2R_arr_x + wq(iq)*dt*(RL2iq_arr_x); %quadrature formula
                    
                    if (L2L2R<0)
                        disp('L2L2R<0')
                    end
                end
            end
            
            %             L2L2_arr_cumulative  = L2L2_arr_cumulative + L2L2R_arr;
            
            if it ==0 && mesh_adapt==0
                c_0_coeff_arr_0 = c_0_coeff_arr_old;
                x_0 = x;
                dist_x_pl_0 = dist_x_pl;
                dist_x_min_0 = dist_x_min;
                error_0 = space_int_vector_rec_weno_nonu(x,dist_x_pl,dist_x_min,uold,0*dt,ex,Lx);

            end
            
            % adaptive step where we use estimator to compute refinement
            % before proceeding to the next step
            
            
            if (showplot && mod(reps,fps)==0)
                N_subplots  = 1;
                N_subplots_x = 1;
                N_subplots_y = 3;
                
                subplot(N_subplots_x,N_subplots_y,1)
                plot(x,u,'ro',x,ex(x,t+dt),'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                %                         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(exponent),'$'],'Interpreter','latex','FontSize',f_size)
                ylabel('$u$','Interpreter','latex','FontSize',f_size)
                xlabel('$x$','Interpreter','latex','FontSize',f_size)
                legend('Numerical','Exact','Location','SouthEast')
                pbaspect([1 1 1])
                grid on;
                ylim([-1 1])
                set(gca, 'Fontsize',f_size);
                
                subplot(N_subplots_x,N_subplots_y,2)
                plot(x,RL2iq_arr,'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                ylabel(['$\mathcal{E}\left(t^n; \left[x_j, x_{j+1}\right]\right)$'],'Interpreter','latex','FontSize',f_size)
                xlabel('$x$','Interpreter','latex','FontSize',f_size)
                pbaspect([1 1 1])
                %                 set(gca, 'Fontsize',f_size);
                
                
                subplot(N_subplots_x,N_subplots_y,3)
                plot(x,dist_x_pl,'b.') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                ylim([0, max(dist_x_pl)])
                ylabel(['$Grid\, Spacing$'],'Interpreter','latex','FontSize',f_size)
                xlabel('$x$','Interpreter','latex','FontSize',f_size)
                %                 set(gca, 'Fontsize',f_size);
                pbaspect([1 1 1])
                
                
                %                 subplot(N_subplots_x,N_subplots_y,4)
                %                 plot(x,dist_x_pl,'b.') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                %                 ylim([0, max(dist_x_pl)])
                %                 ylabel(['$N\left(dofs\right)$'],'Interpreter','latex','FontSize',f_size)
                %                 xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
                %                 set(gca, 'Fontsize',f_size);
                
                
                
                
                
                
                
                pause(0.01)
                set(gcf, 'Position',  [100, 100, 1200, 400]) % position position width height
                
                
                %                 pause(0.01)
                %                 set(gcf, 'Position',  [100, 100, 1000, 300]) % position position width height
                %
                %                 if (reps==20|| reps==3500)
                %                     saveas(figgy_parasite,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/fig_',scheme_arr{i_scheme},'_',IC_string,'_IC_','t=',num2str(t),'_no_mesh_adapt_1D'],'png')
                %                     saveas(figgy_parasite,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/fig_',scheme_arr{i_scheme},'_',num2str(reps),'_no_mesh_adapt_1D_',num2str(i_interpolant)],'png')
                %                     reps
                %                     L2L2R
                %                 end
                
                if (record_video)
                    F_var(rep_video) = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                    
                end
                rep_video = rep_video+1;
                
                
            end
            
            
            
            if it ==0
                % t=0;
                bound_arr(1) = sqrt(exp(0*dt)*error_0);
%                 bound_arr(1) = sqrt(exp(it*dt)*(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
%                 error_new = sqrt(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,u,0*dt,ex,c_0_coeff_arr_0));
%                 error_arr(1) =error_new;% max(error_new, error_old);
%                 error_old = max(error_arr);
                error_new = sqrt(error_0);
                error_arr(1) = error_new;
                error_old = max(error_arr);
                
                EI_index(1) = bound_arr(1)./error_arr(1);
                error_new_arr(end+1) = error_new;
                dofs_arr(1) = length(x_0);
                
                % t=dt
                bound_arr(end+1) = sqrt(exp(1*dt)*(L2L2R + error_0));
%                 bound_arr(end+1) = sqrt(exp(1*dt)*(L2L2R + space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
                error_new = sqrt( space_int_vector_rec_weno_nonu(x,dist_x_pl,dist_x_min,u,1*dt,ex,Lx));
%                 error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,1*dt,ex,c_0_coeff_arr_new));
                error_arr(end+1) =error_new;% max(error_new, error_old);
                error_old = max(error_arr);
                EI_index(end+1) = bound_arr(end)./error_arr(end);
                error_new_arr(end+1) = error_new;
                dofs_arr(end+1) = length(x);
                
                
            else
                % t>dt
                bound_arr(end+1) = sqrt(exp((it+1)*dt)*(L2L2R + error_0));
%                 bound_arr(end+1) = sqrt(exp((it+1)*dt)*(L2L2R + space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
                error_new = sqrt( space_int_vector_rec_weno_nonu(x,dist_x_pl,dist_x_min,u,(it+1)*dt,ex,Lx));
%                 error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,(it+1)*dt,ex,c_0_coeff_arr_new));
                error_arr(end+1) = error_new;%max(error_new, error_old);
                error_old = max(error_arr);
                EI_index(end+1) = bound_arr(end)./error_arr(end);
                error_new_arr(end+1) = error_new;
                dofs_arr(end+1) = length(x);
                
            end
            
            
            
            
            
            
            
            it = it+1;
            t = t + dt; %move in time
            %             %% refine the mesh
            %             M.asgnVals(u, RL2iq_arr);
            %             M.mark_ref_coar();
            %             M.ref_coarsen_local();
            %             [x_ref, u_ref] = M.getNewVals();
            %
            %             x_ref = [x_ref, Lx];
            %             dist_x_pl_ref = x_ref(2:end) -x_ref(1:end-1);
            %             dist_x_min_ref = circshift(dist_x_pl_ref, 1);
            %             dist_x_pl = dist_x_pl_ref;
            %             dist_x_min = dist_x_min_ref;
            %             x = x_ref(1:end-1);
            %
            %
            %             dt = .1 *min(dist_x_pl);
            
            uold = u;
            reps = reps+1;
            time_arr(end+1) = t;
            
        end
        if (record_video)
            
            video = VideoWriter(['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/videos/WENO3_nonu_adapt_linadvect_step_ic','_',num2str(maxit_arr(m)),'_solution_with_',adaptivity,'.mp4'],'MPEG-4');
            video.FrameRate = fps;
            open(video);
            writeVideo(video,F_var);
            close(video);
        end
        
        %
        finalL2err(m) = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,T,ex,c_0_coeff_arr_new)); %compute final time L2 error
        R(m) = sqrt(exp(T)*(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0) + L2L2R)); %bound
        
        EOCe(1) = 0;
        EOCR(1) = 0;
        
        DOFs(m) = length(x);
        %         DOFs(m) = length(x);
        if m > 1
            %             EOCe(m) = log(finalL2err(m-1)/finalL2err(m))/log(2);
            %             EOCR(m) = log(R(m-1)/R(m))/log(2);
            EOCe(m) = log(finalL2err(m-1)/finalL2err(m))/log(DOFs(m)/DOFs(m-1));
            EOCR(m) = log(R(m-1)/R(m))/log(2);
            EOC_error_arr_output(m) = EOCe(m);
            EOC_bound_arr_output(m) = EOCR(m);
        end
        EI(m) = R(m)/finalL2err(m);
        fprintf(1,'||(u - IU)(T)||_L2 = %.5f  EOC = %1.2f\n',finalL2err(m),EOCe(m))
        fprintf(1,'      ||R||_L2(L2) = %.5f  EOC = %1.2f\n',R(m),EOCR(m))
        fprintf(1,'                EI = %.5f  \n',EI(m))
        
        % Arrays for plotting
        
        
        cell_cell_arr{m}=[time_arr;bound_arr;error_arr;EI_index;dofs_arr];
        
        disp(['Iter', num2str(m),'has passed'])
    end
%     save([scheme_arr{i_scheme},'_cell_arr_file_disc_adv_for_thesis_',init_conds,'_',interpolant,'_lin_advect_comparison_',adaptivity,'.mat'],'cell_cell_arr')
    
end


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