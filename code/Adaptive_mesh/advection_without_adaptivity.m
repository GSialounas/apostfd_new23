clear all;
close all;
mesh_adapt=1;
maxit_arr = [1:4]; % max refinement iterations
record_video_arr = zeros(size(maxit_arr));
showplot_arr = zeros(size(maxit_arr));

record_video_arr(end)= 0;
showplot_arr(end) =0;
% record_video = 0;
fps= 50;

rep_video = 1;
f_size = 16;
max_number_ref_levels = 1;
gamma_param = .5;

scheme_arr = {'FTBS'};
intrpl_arr = {'P1'};
init_conds = 'stepIC';
% interpolant_arr
exponent_arr = [1];
global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();
showplot=1;
for i_scheme = 1:length(scheme_arr)
    
    
    cell_cell_arr = {};
    figgy_parasite = figure();
    

    scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
    interpolant = intrpl_arr{1};
    exponent =exponent_arr(1);
    IC_string = 'sin';
    
    T = 1.5; % final time
   
    Lx = 1;
    ratio_cTf = 1;
    part_fine = 1;
    
    ex = @(x,t) sin(2*pi*(x-t)); %smooth exact solution
    ex_x =  @(x,t) 2*pi*cos(2*pi*(x-t));
    
    
    cntr = .25;
    Lx=1;
    pow_exp=100;
    %             ratio_cTf=2;
    part_fine=1;
    radius = 0.125;
    % ex = @(x,t) exp(-100*(x-cntr - t).^2); %smooth exact solution
    ex = @(x,t) fn_stp_exact(x,cntr,t,Lx,radius) ;%
    %  ex = @(x,t) exp(-10*((10*(x-cntr - t))).^2); %smooth exact solution
    %  grid_x_sin=@(x) (sin(50*pi*x)).^2;
    %   ex = @(x,t) max(exp(-pow_exp*((mod(x-cntr-t,Lx)).^2)),exp(-pow_exp*((mod(x-cntr-t,-Lx)).^2))); %smooth exact solution
    
    
    fprintf(1,'Using the %s scheme\n',scheme)
    maxit = length(maxit_arr);

    
    for m = 1:length(maxit_arr)
        record_video = record_video_arr(m);
        %                 showplot =showplot_arr(m);
        
        
        
        h = 2^(-(maxit_arr(m)+6)); % mesh size
        x = create_grid(Lx, ratio_cTf, h, part_fine);
        dt = .1*h^1;
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
        bound_arr = [];%zeros(1,ceil(T/dt));
        error_arr = [];% zeros(1,ceil(T/dt));
        time_arr = [];% zeros(1,ceil(T/dt));
        dofs_arr = [];% zeros(1,ceil(T/dt));
        EI_index = [];
        
        time_arr(1) = 0;
        
        u= uold;
        
        
        L2L2_arr_cumulative= zeros(size(x));
        error_old = 0;
        error_new = 0;
        if m ==1
            dt_coarsest=dt;
        end
        error_new_arr=[0];
        
        while t <= T-dt/2
            
            % initial solve at time step  t = it*dt
            spatial_disc = "BS";
%             if it==0
%                 uold_x = ex_x(x,0);
%             else
                uold_x = 1./(dist_x_min).*(-circshift(uold,1) + circshift(uold,0)); 
%             end
            
            u= uold - dt*uold_x;
            u_x = 1./(dist_x_min).*(-circshift(u,1) + circshift(u,0));

            
            % Initial indicator calculation at time step t = it*dt
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
            L2L2_arr_cumulative  = L2L2_arr_cumulative + L2L2R_arr;
            
            if it ==0
                c_0_coeff_arr_0 = c_0_coeff_arr_old;
                x_0 = x;
                dist_x_pl_0 = dist_x_pl;
                dist_x_min_0 = dist_x_min;
            end
            
            
            if (showplot && mod(reps,fps)==0)
                N_subplots  = 1;
                N_subplots_x = 1;
                N_subplots_y = 2;
                
                subplot(N_subplots_x,N_subplots_y,1)
                plot(x,u,'ro',x,ex(x,t),'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                %                         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(exponent),'$'],'Interpreter','latex','FontSize',f_size)
                ylabel('$u$','Interpreter','latex','FontSize',f_size)
                legend('Numerical','Exact','Location','SouthEast')
                pbaspect([1 1 1])
                grid on;
                ylim([-1 1])
                set(gca, 'Fontsize',f_size);
                
                subplot(N_subplots_x,N_subplots_y,2)
                plot(x,RL2iq_arr,'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                ylabel(['$\mathcal{E}\left(t^n; \left[x_j, x_{j+1}\right]\right)$'],'Interpreter','latex','FontSize',f_size)

                pbaspect([1 1 1])

                
                pause(0.01)
                set(gcf, 'Position',  [100, 100, 1000, 500]) % position position width height
                
                
                %                         pause(0.01)
                %                         set(gcf, 'Position',  [100, 100, 1000, 300]) % position position width height
                %
                %                         if (reps==20|| reps==3500)
                %                             saveas(figgy_parasite,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/fig_',scheme_arr{i_scheme},'_',IC_string,'_IC_','t=',num2str(t),'_no_mesh_adapt_1D'],'png')
                %                             saveas(figgy_parasite,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/fig_',scheme_arr{i_scheme},'_',num2str(reps),'_no_mesh_adapt_1D_',num2str(i_interpolant)],'png')
                %                             reps
                %                             L2L2R
                %                         end
                
                if (record_video)
                    F_var(rep_video) = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                    
                end
                rep_video = rep_video+1;
                
                
            end
            
            
            %                     bound_arr_eoc(it) = sqrt(L2L2R*exp(it*dt) + space_int_vector(x,dist_x_pl,dist_x_min,ex(x,0),0,ex)); % the second part of the expression is the initial error(e_0)
            
            if it ==0
                % t=0;
                bound_arr(1) = sqrt(exp(it*dt)*(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
                error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,0*dt,ex,c_0_coeff_arr_0));
                error_arr(1) = max(error_new, error_old);
                error_old = max(error_arr);
                EI_index(1) = bound_arr(1)./error_arr(1);
                error_new_arr(end+1) = error_new;
                
                
                % t=dt
                bound_arr(end+1) = sqrt(exp(1*dt)*(L2L2R + space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
                error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,1*dt,ex,c_0_coeff_arr_new));
                error_arr(end+1) = max(error_new, error_old);
                error_old = max(error_arr);
                EI_index(end+1) = bound_arr(end)./error_arr(end);
                error_new_arr(end+1) = error_new;
                
                
            else
                % t>dt
                bound_arr(end+1) = sqrt(exp((it+1)*dt)*(L2L2R + space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
                error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,(it+1)*dt,ex,c_0_coeff_arr_new));
                error_arr(end+1) = max(error_new, error_old);
                error_old = max(error_arr);
                EI_index(end+1) = bound_arr(end)./error_arr(end);
                error_new_arr(end+1) = error_new;
                
            end
            
            
            it = it+1;
            t = t + dt; %move in time
            
            uold = u;
            reps = reps+1;
            time_arr(end+1) = t;
            
        end
        if (record_video)
            video = VideoWriter(['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/FTBS_uniform_step_ic','_',num2str(maxit_arr(m)),'_solution.mp4'],'MPEG-4');
            video.FrameRate = fps;
            open(video);
            writeVideo(video,F_var);
            close(video);
        end
        
        
        finalL2err(m) = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,T,ex,c_0_coeff_arr_new)); %compute final time L2 error
        R(m) = sqrt(exp(T)*(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0) + L2L2R)); %bound
        
        EOCe(1) = 0;
        EOCR(1) = 0;
        if m > 1
            EOCe(m) = log(finalL2err(m-1)/finalL2err(m))/log(2);
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
        
    end
%     save([scheme_arr{i_scheme},'_cell_arr_file_',init_conds,'_',interpolant,'_lin_advect.mat'],'cell_cell_arr')
    
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