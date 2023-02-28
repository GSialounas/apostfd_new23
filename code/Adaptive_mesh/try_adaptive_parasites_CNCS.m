clear all;
close all;
clc;

record_video = 0;
showplot=1;

mesh_adapt=0;
if mesh_adapt == 1
    adaptivity = 'adapt_ON';
else
    adaptivity = 'adapt_OFF';
end
maxit_arr = [1:4]; % max refinement iterations

fps= 20;

rep_video = 1;
f_size = 16;
max_number_ref_levels = 1;
gamma_param = .5;

% scheme_arr = {'FTBS'};
scheme_arr = {'CNCS'};

intrpl_arr = {'P2'};
init_conds = 'expIC';
% interpolant_arr
exponent_arr = [1];
global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();
for i_scheme = 1:length(scheme_arr)
    
    
    cell_cell_arr = {};
    figgy_parasite = figure();
    
    
    scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
    interpolant = intrpl_arr{1};
    exponent =exponent_arr(1);
    IC_string = 'sin';
    
    T = .5; % final time
    U=1;
    K=0;
    theta = .5;
    Lx = 1;
    ratio_cTf = 2;
    part_fine = .5;
    
    %     ex = @(x,t) sin(2*pi*(x-t)); %smooth exact solution
    %     ex_x =  @(x,t) 2*pi*cos(2*pi*(x-t));
    
    
    cntr = .25;
    Lx=1;
    radius = 0.125;
    % ex = @(x,t) exp(-100*(x-cntr - t).^2); %smooth exact solution
%     ex = @(x,t) fn_stp_exact(x,cntr,t,Lx,radius) ;%
    ex = @(x,t) exp(-10*((10*(x-cntr - t))).^2);
    %  ex = @(x,t) exp(-10*((10*(x-cntr - t))).^2); %smooth exact solution
    %  grid_x_sin=@(x) (sin(50*pi*x)).^2;
    %   ex = @(x,t) max(exp(-pow_exp*((mod(x-cntr-t,Lx)).^2)),exp(-pow_exp*((mod(x-cntr-t,-Lx)).^2))); %smooth exact solution
    
    
    fprintf(1,'Using the %s scheme\n',scheme)
    maxit = length(maxit_arr);
    
    
    for m =length(maxit_arr):length(maxit_arr)
        %                 showplot =showplot_arr(m);
        
        
        % produce initial mesh
        h = 2^(-(maxit_arr(m)+6)); % mesh size
        %         dt_coarsest = .1 * 2^(-(maxit_arr(end)+8)); % mesh size
        if mesh_adapt ==0
%             x= linspace(0,1,741);
            x = create_grid(Lx, ratio_cTf, h, part_fine);
            
        else
            x = create_grid(Lx, ratio_cTf, h, part_fine);
            
        end
        dt =  .1 *h;%.1*h^1; % initial time-step
        dx_arr(m) = h;
        dt_arr(m) = dt; % timestep size
        
        x=x(1:end-1);
        dist_x_pl = circshift(x,-1)-x;
        dist_x_pl(end)= Lx-x(end);
        dist_x_min = x- circshift(x,1);
        dist_x_min(1)= Lx-x(end);
        
        
        
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
        
        % form initial mesh
        tree_list = [];
        root_ids = [];
        x = [x, Lx];
        tol = 1e-6;
        max_ref= 4;
        for i = 1:length(x)-1
            if i==1
                tree_list = [ Tree(element(x(i), x(i+1),i,tol),i,tol)];
                %                 tree_list(i) = tree_list(i).assign_vals([sin(x(i)), sin(x(i+1))], exp(-100*(x(i)-.5).^2));
                root_ids(i) = i;
            else
                tree_list(end+1) = Tree(element(x(i), x(i+1),i,tol),i,tol);
                %                 tree_list(i) = tree_list(i).assign_vals([sin(x(i)), sin(x(i+1))], exp(-100*(x(i)-.5).^2));
                root_ids(i) = i;
            end
        end
        
        if mesh_adapt==1
            M = mesh(tree_list,root_ids,max_ref);
            M.updateMesh();
            x = x(1:end-1);
            M.asgnVals(uold, zeros(size(uold)));
            
            M.refineGlobal();
            M.refineGlobal();
            M.refineGlobal();
            M.refineGlobal();
            [x_ref,uold_ref] = M.getNewVals();
            
            % amend the grid and uld
            x_ref = [x_ref, Lx];
            uold  = uold_ref;
            
            dist_x_pl_ref = x_ref(2:end) -x_ref(1:end-1);
            dist_x_min_ref = circshift(dist_x_pl_ref, 1);
            dist_x_pl = dist_x_pl_ref;
            dist_x_min = dist_x_min_ref;
            x = x_ref(1:end-1);
        else
            x=x(1:end-1);
        end
        if mesh_adapt == 1
            uold= ex(x,0);
            M.asgnVals(uold, zeros(size(uold)));
        end
        
        
        stp_chgmat = 0;
%  
        A_adv  = createA_advection_TRISTAN_nonu_central(U,dt, x,Lx);%createA_advection(U,dt, x);
        A_diff = createA_diffusion_TRISTAN(dt,x, K,Lx);%createA_diffusion(dt,x, K);
        I   = speye(size(A_adv));
        Nt  = (floor(T/dt)); % Number of time-steps

        shift_mat = ( I + theta * (A_adv - stp_chgmat*A_diff) )\( (I + (1-theta) * (-A_adv + stp_chgmat*A_diff) ) );
        shift_mat = shift_mat^1;
%         
        
        while t <= T-dt/2          
            % initial solve at time step  t = it*dt
            spatial_disc = "CS";
            
            % Initial solve and initial estimator computation
            uold_x = fn_central_nonu(dist_x_pl, dist_x_min, uold);
%             uold_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(uold,1) + circshift(uold,-1));
%             uold = uold';
%             u = uold -dt*uold_x;%shift_mat * uold;
            u = shift_mat * uold';

            u = u';
%             uold = uold';            
%             u_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(u,1) + circshift(u,-1));
            u_x = fn_central_nonu(dist_x_pl, dist_x_min, uold);
         
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
                    else
                    end
                    L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                    L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
                    L2L2R_arr_t = RL2iq_arr_t;%L2L2R_arr_t + wq(iq)*dt*(RL2iq_arr_t); %quadrature formula
                    L2L2R_arr_x = RL2iq_arr_x;% L2L2R_arr_x + wq(iq)*dt*(RL2iq_arr_x); %quadrature formula
                    
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
                    elseif intrpl_arr{i_interpolant}=='P2'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_2(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    else
                    end
%                     L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                    L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
                    L2L2R_arr_t = RL2iq_arr_t;%L2L2R_arr_t + wq(iq)*dt*(RL2iq_arr_t); %quadrature formula
                    L2L2R_arr_x = RL2iq_arr_x;% L2L2R_arr_x + wq(iq)*dt*(RL2iq_arr_x); %quadrature formula
                    
                    if (L2L2R<0)
                        disp('L2L2R<0')
                    end
                end
                if it ==0 && mesh_adapt==1
                    c_0_coeff_arr_0 = c_0_coeff_arr_old;
                    x_0 = x;
                    dist_x_pl_0 = dist_x_pl;
                    dist_x_min_0 = dist_x_min;
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
   
        
%                 A_adv  = createA_advection_TRISTAN_nonu_central(U,dt, x,Lx);%createA_advection(U,dt, x);
%                 A_diff = createA_diffusion_TRISTAN(dt,x, K,Lx);%createA_diffusion(dt,x, K);
%                 I   = speye(size(A_adv));
%                 Nt  = (floor(T/dt)); % Number of time-steps
%                 
%                 shift_mat = ( I + theta * (A_adv - stp_chgmat*A_diff) )\( (I + (1-theta) * (-A_adv + stp_chgmat*A_diff) ) );
%                 shift_mat = shift_mat;
                 uold_x = fn_central_nonu(dist_x_pl, dist_x_min, uold);
%                 uold_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(uold,1) + circshift(uold,-1));
%                 uold = uold';
%                 u = shift_mat * uold;
                u = uold -dt*uold_x;
%                 u = u';
%                 uold = uold';
%                 u_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(u,1) + circshift(u,-1));
                 u_x = fn_central_nonu(dist_x_pl, dist_x_min, u);
                % bound for plotting purposes
                
                L2L2R_arr= zeros(size(x));
                L2L2R_arr_t= zeros(size(x));
                L2L2R_arr_x= zeros(size(x));
                for iq = 1 : nq
                    tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                    if interpolant=='P1'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_1(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    elseif intrpl_arr{i_interpolant}=='P2'
                        [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,RL2iq_arr_t,RL2iq_arr_x] = compute_Rs_vector_temp_1_spatiotemp_2(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, uold_x, u_x,spatial_disc); %compute R at gauss point
                    else
                    end
                    L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                    L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
                    L2L2R_arr_t = RL2iq_arr_t;%L2L2R_arr_t + wq(iq)*dt*(RL2iq_arr_t); %quadrature formula
                    L2L2R_arr_x = RL2iq_arr_x;% L2L2R_arr_x + wq(iq)*dt*(RL2iq_arr_x); %quadrature formula
                    
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
            end
            
            % adaptive step where we use estimator to compute refinement
            % before proceeding to the next step
            
            
            if (showplot && mod(reps,fps)==0)
                N_subplots  = 1;
                N_subplots_x = 1;
                N_subplots_y = 3;
                
                subplot(N_subplots_x,N_subplots_y,1)
%                 plot(x,u,'ro',x,ex(x,t+dt),'b')
                plot(x,u,'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')

                %                         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(exponent),'$'],'Interpreter','latex','FontSize',f_size)
                ylabel('$u$','Interpreter','latex','FontSize',f_size)
                xlabel('$x$','Interpreter','latex','FontSize',f_size)
%                 legend('Numerical','Exact','Location','SouthEast')
                pbaspect([1 1 1])
                grid on;
                ylim([-2e-6 2e-6])
                set(gca, 'Fontsize',f_size);
                
                subplot(N_subplots_x,N_subplots_y,2)
                plot(x,RL2iq_arr,'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                ylabel(['$\mathcal{E}\left(t^n; \left[x_j, x_{j+1}\right]\right)$'],'Interpreter','latex','FontSize',f_size)
                xlabel('$x$','Interpreter','latex','FontSize',f_size)
                pbaspect([1 1 1])
                ylim([-2e-6 2e-6])

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
                bound_arr(1) = sqrt(exp(it*dt)*(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
                error_new = sqrt(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,u,0*dt,ex,c_0_coeff_arr_0));
                error_arr(1) =error_new;% max(error_new, error_old);
                error_old = max(error_arr);
                EI_index(1) = bound_arr(1)./error_arr(1);
                error_new_arr(end+1) = error_new;
                dofs_arr(1) = length(x_0);
                
                % t=dt
                bound_arr(end+1) = sqrt(exp(1*dt)*(L2L2R + space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
                error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,1*dt,ex,c_0_coeff_arr_new));
                error_arr(end+1) =error_new;% max(error_new, error_old);
                error_old = max(error_arr);
                EI_index(end+1) = bound_arr(end)./error_arr(end);
                error_new_arr(end+1) = error_new;
                dofs_arr(end+1) = length(x);
                
                
            else
                % t>dt
                bound_arr(end+1) = sqrt(exp((it+1)*dt)*(L2L2R + space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
                error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,(it+1)*dt,ex,c_0_coeff_arr_new));
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
            
            video = VideoWriter(['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/videos/CNCS_para_adapt_linadvect_exp_ic','_',num2str(maxit_arr(m)),'_solution_with_',adaptivity,'.mp4'],'MPEG-4');
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
%     save([scheme_arr{i_scheme},'_cell_arr_file_',init_conds,'_',interpolant,'_lin_advect_comparison_',adaptivity,'.mat'],'cell_cell_arr')
    
end


function [L2Rt,L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old,L2Rt_arr_t,L2Rt_arr_x] = compute_Rs_vector_temp_1_spatiotemp_2(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
L2Rt_arr = zeros(size(x));
L2Rt_arr_t = zeros(size(x));
L2Rt_arr_x = zeros(size(x));
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients

% reconstruction with derivative matching at t_n
c_0_t = uold; % this is u at t_n
c_1_t = (1/dt)*(u-uold);

diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t * diff_t ;
Ut_t = c_1_t;

% We calculate the spatial derivative discretisation using the spatial
% reconstruction at time t
dx_fine = x(2)-x(1);
if spatial_disc == 'CS'
    Ut_x = fn_central_nonu(dist_x_pl, dist_x_min, Ut);
    Ut_xt = fn_central_nonu(dist_x_pl, dist_x_min, Ut_t);
    %         Ut_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut,1) + circshift(Ut,-1));
    %         Ut_xt = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut_t,1) + circshift(Ut_t,-1));
elseif spatial_disc == 'BS'
    Ut_x  = 1./dist_x_pl.*(circshift(Ut,-1)-circshift(Ut,0));%1./(dist_x_pl).*(flux_lxf(dist_x_pl, Ut, dt)-flux_lxf(dist_x_min,circshift(Ut,1),dt));% fn_backward_burger(dist_x_pl, dist_x_min, Ut);
    Ut_xt = 1./dist_x_pl.*(circshift(Ut_t,-1)-circshift(Ut_t,0));%1./dist_x_pl.*(flux_lxf(dist_x_pl, Ut_t, dt)-flux_lxf(dist_x_min,circshift(Ut_t,1),dt));  %fn_backward_burger(dist_x_pl, dist_x_min, Ut_t);
elseif spatial_disc == '2S'
    Ut_x  = fn_2S(dist_x_pl,Ut);
    Ut_xt = fn_2S(dist_x_pl,Ut_t);
    %     Ut_x= 1./(2*h)*( 3*circshift(Ut,0)  - 4*circshift(Ut,1) + circshift(Ut,2) );
    %     Ut_xt= 1./(2*h)*( 3*circshift(Ut_t,0)  - 4*circshift(Ut_t,1) + circshift(Ut_t,2) );
elseif spatial_disc == '3S'
    Ut_x = 1./(6*h)*( 2*circshift(Ut,-1)  + 3*circshift(Ut,0) - 6* circshift(Ut,1) +circshift(Ut,2));
    Ut_xt = 1./(6*h)*( 2*circshift(Ut_t,-1)  + 3*circshift(Ut_t,0) - 6* circshift(Ut_t,1) +circshift(Ut_t,2));
elseif spatial_disc == '4CS'
    Ut_x = -1./(12*h).*(-circshift(Ut,2)+8*circshift(Ut,1) - 8*circshift(Ut,-1) + circshift(Ut,-2));
    Ut_xt = -1./(12*h).*(-circshift(Ut_t,2)+8*circshift(Ut_t,1) - 8*circshift(Ut_t,-1) + circshift(Ut_t,-2));
else
end

% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
% Reconstruction where we match the spatial derivative at x_j
c_0_ts = Ut; % this is u at t_n
c_1_ts = Ut_x; % - f_h(U^n)
c_2_ts = (1./(dist_x_pl.^2)) .* (circshift(Ut,-1) - (Ut + dist_x_pl .* Ut_x));

c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = Ut_xt; % - f_h(U^n)
c_2_ts_t = (1./(dist_x_pl.^2)) .* (circshift(Ut_t,-1) - (Ut_t + dist_x_pl .* Ut_xt));



% c_old and c_new with matching derivative at x_j
c_0_old = uold;
c_1_old = uold_x;
c_2_old = (1./(dist_x_pl.^2)) .* (circshift(uold,-1) - (uold + dist_x_pl.* uold_x));

c_0_new = u;
c_1_new = u_x;
c_2_new = (1./(dist_x_pl.^2)) .* (circshift(u,-1) - (u + dist_x_pl.*u_x));

c_0_coeff_arr_old  =  [c_0_old; c_1_old; c_2_old];
c_0_coeff_arr_new  =  [c_0_new; c_1_new; c_2_new];


for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU = c_0_ts + c_1_ts .* diff_x + c_2_ts .* diff_x .^2 ;
    IUx = ( c_1_ts + 2* c_2_ts .*diff_x);
    IUt = (c_0_ts_t + c_1_ts_t.*diff_x +  c_2_ts_t.*diff_x.^2 );
    
    integral_txq = (wq(iq)*dist_x_pl.*( (IUt+IUx)).^2);
    integral_txq_t = (wq(iq)*dist_x_pl.*( IUt).^2);
    integral_txq_x = (wq(iq)*dist_x_pl.*( IUx).^2);
    
    L2Rt = L2Rt + sum(integral_txq);
    L2Rt_arr = L2Rt_arr + (integral_txq);
    L2Rt_arr_t =  c_1_ts;%c_2_ts_t.*diff_x.^2 ;% L2Rt_arr_t + (integral_txq_t);
    L2Rt_arr_x =  c_1_ts_t; %c_2_ts_t.*diff_x.^2 ;%L2Rt_arr_x + (integral_txq_x);
    
    
end
% end

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