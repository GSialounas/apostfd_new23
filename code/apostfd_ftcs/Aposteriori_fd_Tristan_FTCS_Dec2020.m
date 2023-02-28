clear all;
close all;

mesh_adapt=0;
maxit_arr = [1:3]; % max refinement iterations
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

scheme_arr = {'FTCS'};
intrpl_arr = {'P2'};
% interpolant_arr
exponent_arr = [1];
global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();
showplot=0;
for i_scheme = 1:length(scheme_arr)
    
    
    cell_cell_arr = {};
    figgy_parasite = figure();
    
    for i_exponent = 1:length(exponent_arr)
        
        
        for i_interpolant = 1:length(intrpl_arr)
            
            
            scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
            interpolant = intrpl_arr{i_interpolant};
            exponent =exponent_arr(i_exponent);
            IC_string = 'sin';
            
            T = .05; % final time
            U=1;
            K=0;
            
            
            transform_factor = 1;
            cntr = .2;
            Lx=1;
            pow_exp=1000;
            
            %             ratio_cTf_arr = [.5,.5,1,1]
            
            part_fine_arr = [.5,.5,1,1];
            ratio_cTf = 2;
            part_fine = 1;
            
            ex = @(x,t) sin(2*pi*(x-t)); %smooth exact solution
            ex_x =  @(x,t) 2*pi*cos(2*pi*(x-t));
            %           ex = @(x,t) exp(-10*((10*(x-cntr - t))).^2); %smooth exact solution
            %           grid_x_sin=@(x) (sin(50*pi*x)).^2;
            %           ex = @(x,t) max(exp(-pow_exp*((mod(x-cntr-t,Lx)).^2)),exp(-pow_exp*((mod(x-cntr-t,-Lx)).^2))); %smooth exact solution
            
            
            fprintf(1,'Using the %s scheme\n',scheme)
            maxit = length(maxit_arr);
            
            %             ratio_cTf = 2;
            time1= zeros(1,maxit);
            time2 =zeros(1,maxit);
            dx_arr = zeros(1,maxit);
            dt_arr = zeros(1,maxit);
            
            error_arr_output = zeros(1,maxit);
            bound_arr_output = zeros(1,maxit);
            EOC_error_arr_output = zeros(1,maxit);
            EOC_bound_arr_output = zeros(1,maxit);
            
            
            for m = 1:length(maxit_arr)
                record_video = record_video_arr(m);
                %                 showplot =showplot_arr(m);
                
                
                
                h = 2^(-(maxit_arr(m)+3)); % mesh size
                x = create_grid(Lx, ratio_cTf, h, part_fine);
                dt = .2*h^3;%exponent;
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
                
                
                ref_level = zeros(size(x));
                mark_ele_toref = zeros(size(x));
                L2L2_arr_cumulative= zeros(size(x));
                error_old = 0;
                error_new = 0;
                if m ==1
                    dt_coarsest=dt;
                end
                error_new_arr=[0];
                
                while t <= T-dt/2
                    
                    % initial solve at time step  t = it*dt
                    spatial_disc = 'CS';
                    if it==0
                        uold_x = ex_x(x,0);
                    else
                        uold_x = fn_central_nonu(dist_x_pl, dist_x_min, uold);%1./(dist_x_min+dist_x_pl).*(-circshift(uold,1) + circshift(uold,-1));
                    end
                    
                    u= uold - dt*uold_x;
                    u_x = fn_central_nonu(dist_x_pl, dist_x_min, u);%;1./(dist_x_min+dist_x_pl).*(-circshift(u,1) + circshift(u,-1));
                    %                     u_x = 1./(dist_x_min).*(-circshift(u,1) + circshift(u,0));
                    
                    
                    % Initial indicator calculation at time step t = it*dt
                    L2L2R_arr= zeros(size(x));
                    L2L2R_arr_t= zeros(size(x));
                    L2L2R_arr_x= zeros(size(x));
                    for iq = 1 : nq
                        tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                        if intrpl_arr{i_interpolant}=='P1'
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
                    
                    
                    rad_tbl = .1;  pow_sp = 100; h_tbl = 1e-3; tol = 1e-8;
                    stp_chg = zeros(size(x));
                    
                    x_cntrs = x( L2L2R_arr>tol);
                    if (isempty(x_cntrs)==0)
                        stp_chg = make_tbl_sponge_vector(x, stp_chg, h_tbl, rad_tbl, x_cntrs',pow_sp);
                        
                    end
                    if (showplot && mod(reps,fps)==0)
                        N_subplots  = 1;
                        N_subplots_x = 2;
                        N_subplots_y = 1;
                        %                         subplot(N_subplots_x,N_subplots_y,1)
                        plot(x,u,'ro',x,ex(x,t),'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                        %                         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(exponent),'$'],'Interpreter','latex','FontSize',f_size)
                        ylabel('$u$','Interpreter','latex','FontSize',f_size)
                        legend('Numerical','Exact','Location','SouthEast')
                        pbaspect([1 1 1])
                        grid on;
                        ylim([-1 1])
                        set(gca, 'Fontsize',f_size);
                        pause(0.01)
                        set(gcf, 'Position',  [100, 100, 400, 400]) % position position width height
                        
                        
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
                        error_arr(end+1) = error_new;%max(error_new, error_old);
                        error_old = error_new;%max(error_arr);
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
                    video = VideoWriter(['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/FTBS_uniform_sin_ic','_',num2str(maxit_arr(m)),'_solution.mp4'],'MPEG-4');
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
                
                
                cell_cell_arr{i_exponent}{i_interpolant}{m}=[time_arr;bound_arr;error_arr;EI_index;dofs_arr];
                
            end
            
            
            
            
            
        end
        
    end
    %     set(gcf, 'Position',  [100, 100, 1000, 500]) % position position width height
    
    %     saveas(figgy,['fig_',scheme_arr{i_scheme},' L2_err_bound_',IC_string,' ','_Tristan_nonu_middle_parasite'],'png')
    interpolant = intrpl_arr{i_interpolant};
    %     save([scheme_arr{i_scheme},'_cell_arr_file_sin_IC_uniform.mat'],'cell_cell_arr')
    save([scheme_arr{i_scheme},'_cell_arr_file_sin_IC_uniform_',interpolant,'.mat'],'cell_cell_arr')
    
end

function [L2Rt,L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old,L2Rt_arr_t,L2Rt_arr_x] = compute_Rs_vector_temp_1_spatiotemp_1(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc)
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
% h=dist_x_pl(1);
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
c_1_ts = 1./dist_x_pl.*(circshift(Ut,-1)-circshift(Ut,0));%Ut_x; % - f_h(U^n)


c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = 1./dist_x_pl.*(circshift(Ut_t,-1)-circshift(Ut_t,0));%Ut_xt; % - f_h(U^n)


% c_old and c_new with matching derivative at x_j
c_0_old = uold;
c_1_old = 1./dist_x_pl.*(circshift(uold,-1)-uold);%uold_x;

c_0_new = u;
c_1_new = 1./dist_x_pl.*(circshift(u,-1)-u);

c_0_coeff_arr_old  =  [c_0_old; c_1_old];
c_0_coeff_arr_new  =  [c_0_new; c_1_new];


for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU  = c_0_ts + c_1_ts .* diff_x ;
    IUx =  c_1_ts ;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x;
    
    integral_txq = (wq(iq)*dist_x_pl.*( (IUt+IUx)).^2);
    integral_txq_t = (wq(iq)*dist_x_pl.*( IUt).^2);
    integral_txq_x = (wq(iq)*dist_x_pl.r*( IUx).^2);
    
    L2Rt = L2Rt + sum(integral_txq);
    L2Rt_arr = L2Rt_arr + (integral_txq);
    L2Rt_arr_t =  c_1_ts;%c_2_ts_t.*diff_x.^2 ;% L2Rt_arr_t + (integral_txq_t);
    L2Rt_arr_x =  c_1_ts_t; %c_2_ts_t.*diff_x.^2 ;%L2Rt_arr_x + (integral_txq_x);
    
    
end
% end

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


function [L2Rt,L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old,L2Rt_arr_t,L2Rt_arr_x] = compute_Rs_vector_temp_1_spatiotemp_3(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc)
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
% h=dist_x_pl(1);
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
c_2_ts = 3./(dist_x_pl.^2) .* (circshift(Ut,-1) - Ut) - (1./dist_x_pl) .* (circshift(Ut_x,-1) + 2 * Ut_x);
c_3_ts = 1./(dist_x_pl.^2) .* (circshift(Ut_x,-1) + Ut_x) - 2./(dist_x_pl.^3) .* (circshift(Ut,-1) - Ut);

c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = Ut_xt; % - f_h(U^n)
c_2_ts_t = 3./(dist_x_pl.^2) .* (circshift(Ut_t,-1) - Ut_t) - (1./dist_x_pl) .* (circshift(Ut_xt,-1) + 2 * Ut_xt);
c_3_ts_t = 1./(dist_x_pl.^2) .* (circshift(Ut_xt,-1) + Ut_xt) - 2./(dist_x_pl.^3) .* (circshift(Ut_t,-1) - Ut_t);



c_0_old = uold;
c_1_old = uold_x;
c_2_old = 3./(dist_x_pl.^2) .* (circshift(uold,-1) - uold) - (1./dist_x_pl).*(circshift(uold_x,-1) + 2 * uold_x);
c_3_old = 1./(dist_x_pl.^2) .* (circshift(uold_x,-1) + uold_x) - 2./(dist_x_pl.^3) .* (circshift(uold,-1) - uold);


c_0_new = u;
c_1_new = u_x;
c_2_new = 3./(dist_x_pl.^2) .* (circshift(u,-1)-u) - (1./dist_x_pl).*(circshift(u_x,-1) + 2 * u_x);
c_3_new = 1./(dist_x_pl.^2) .* (circshift(u_x,-1) + u_x) - (2./dist_x_pl.^3) .* (circshift(u,-1) - u);

c_0_coeff_arr_old= [c_0_old; c_1_old; c_2_old; c_3_old];
c_0_coeff_arr_new=  [c_0_new; c_1_new; c_2_new; c_3_new];

%     xl= zeros(size(nq));
for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU = c_0_ts + c_1_ts .* diff_x + c_2_ts .* diff_x .^2 + c_3_ts .* diff_x.^3;
    IUx = c_1_ts + 2* c_2_ts .*diff_x + 3* c_3_ts .*diff_x.^2;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x +  c_2_ts_t.*diff_x.^2 + c_3_ts_t.*diff_x.^3;
    
    
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

function int = space_int_vector_gs(x,dist_x_pl,dist_x_min,u,evalt,ex,c_0_coeff_arr) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq
% c_coefficients = [c_0; c_1; c_2; c_3];
% h = x(2)-x(1);
int = 0;
% for j = 1 : length(x)-1
l_coeffs = length(c_0_coeff_arr);

% diff_x_arr = zeros(size(c_0_coeff_arr));
for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
    diff_x = xiq-x;
    
    IU = sum(c_0_coeff_arr.*(diff_x.^([0:(size(c_0_coeff_arr,1)-1)]')),1);
    
    int = int + sum(wq(iq).* dist_x_pl .* (IU - ex(xiq,evalt)).^2);
end
% end

end

function int = space_int_vector_for_plot(x,dist_x_pl,dist_x_min,u,evalt,ex,c_0_coeff_arr) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq
% c_coefficients = [c_0; c_1; c_2; c_3];
% h = x(2)-x(1);
int = zeros(size(x));
% for j = 1 : length(x)-1
l_coeffs = length(c_0_coeff_arr);

% diff_x_arr = zeros(size(c_0_coeff_arr));

for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
    diff_x = xiq-x;
    %     for ix = 1:size(c_0_coeff_arr,1)
    %         diff_x_arr(ix,:) = diff_x.^(ix-1);
    %     end
    IU = sum(c_0_coeff_arr.*[ones(1,l_coeffs); diff_x; diff_x.^2],1);
    
    
    int = int + (wq(iq).* dist_x_pl .* (IU - ex(xiq,evalt)).^2);
end
% end

end



function mark_ele_toref = mark_ref_coars_max(arr_indic,gamma_param)
% This function marks for refinment based on the maximum strategy (see
% Alberta20 notes pg.32)
% ========== DESCRIPTION ==========
% This function implements the maximum refinement strategy from Alberta20
% pg. 32, whereby the refinement criterion is whether or not the indicator
% for a given element exceeds the maximum value of the indicator on the
% current triangulation times a parameter gamma.
%
% ========== INPUT ==========
% arr_indic   : array containing local indicators (i.e. on each element)
% gamma_param : the refinement critierion as a percentage of the maximum
%               error indicator (see Alberta20 pg. 32)
% ========== OUTPUT ==========
% mark_ele_toref : array who has the same number of elements as the current
%                  triangulation BEFORE ANY REFINEMENT/COARSENING.  If the
%                  index in a given position is 0 - do nothing.  1- refine
%                  and -1 coarsen (coarsening is not implemented yet).
% ========== START OF FUNCTION ==========

mark_ele_toref = zeros(size(arr_indic));
max_ind = gamma_param * max(arr_indic);
for i  = 1 :length(arr_indic)
    if arr_indic(i)>max_ind
        mark_ele_toref(i)= 1;
    end
end
end


function mark_ele_toref = mark_ref_coars_equi(arr_indic,theta_param,tol)

% ========== DESCRIPTION ==========
% This function implements the equid. refinement strategy from Alberta20
% pg. 32, whereby the refinement criterion is whether or not the indicator
% for a given element exceeds theta*tol/N_k^1/p, where N_k is the number of
% elements in the current triangulation.
%
% ========== INPUT ==========
% arr_indic   : array containing local indicators (i.e. on each element)
% theta_param : the refinement critierion as a percentage of the
%               equidistribution error indicator (see Alberta20 pg. 32)
% tol         : chosen tolerance
% ========== OUTPUT ==========
% mark_ele_toref : array who has the same number of elements as the current
%                  triangulation BEFORE ANY REFINEMENT/COARSENING.  If the
%                  index in a given position is 0 - do nothing.  1- refine
%                  and -1 coarsen (coarsening is not implemented yet).
% ========== START OF FUNCTION ==========
N_k = length(arr_indic);
mark_ele_toref = zeros(size(arr_indic));
equid_ind = theta_param * tol/(sqrt(N_k));
for i  = 1 :length(arr_indic)
    if arr_indic(i)>equid_ind
        mark_ele_toref(i) = 1;
    end
end
end

function [x_new, uold_nu_new, uold_x_nu_new, ref_level_new]= fn_prolong(mark_ele_toref,ref_level_old,arr_ele,x,uold_nu, uold_x_nu,Lx,max_number_ref_levels)

% ========== DESCRIPTION ==========
% This function returns the grid after one refinement cycle
% ========== INPUT ==========
% mark_ele_toref  :  vector with same number of entries as elemetns in the
%                    current triangulation.  An entry number of 1 means
%                    refine this element.  An entry number of 0 means do
%                    nothing and -1 means coarsen (not implemented yet)
% arr_ele         :  In the future this will be a binary tree
%                    representation of the mesh.  We might need this to
%                    coarsen in 1D and definetely in 2D.
% Lx              :  This is actually entry x(end), which we don't include
%                    because of periodic boundary conditions but which we
%                    nonetheless need.
% ref_level_old   :  The refinement level of the element compared to its
%                    macro element parent.  A refinement level of 1 means
%                    the element has been refined once from its parent.
% ========== OUTPUT ==========
% x_new           : the new grid after the refinement
% ref_level_new   : the refinement_level after the refinement process is
%                   finished
interp_local  = [1,0;.5,.5;0,1];
x_new=[];
uold_nu_new=[];
uold_x_nu_new=[];
ref_level_new = [];
for i =1:length(mark_ele_toref)
    
    ref_level_ele = ref_level_old(i);
    
    if ((mark_ele_toref(i) ==1)&&ref_level_ele<max_number_ref_levels)
        % refinement produces two new elemets
        if i<length(x)
            ele_new = (interp_local*[x(i);x(i+1)])';
            uold_nu_new_ele = (interp_local*[uold_nu(i);uold_nu(i+1)])';
            uold_x_nu_new_ele = (interp_local*[uold_x_nu(i);uold_x_nu(i+1)])';
        elseif i==length(x)
            ele_new = (interp_local*[x(i);Lx])';
            uold_nu_new_ele = (interp_local*[uold_nu(i);uold_nu(1)])';
            uold_x_nu_new_ele = (interp_local*[uold_x_nu(i);uold_x_nu(1)])';
            % due to periodic boundary conditions)
            ele_new = ele_new(1:end-1);
            uold_nu_new_ele=uold_nu_new_ele(1:end-1);
            uold_x_nu_new_ele = uold_x_nu_new_ele(1:end-1);
            
            %             ele_new = (interp_local*[x(i);Lx])';
            %             % periodic boundary conditions means x_0 and x(end) are
            %             % identified.
            %             uold_nu_new = (interp_local*[uold_nu(i);uold_nu(1)])';
            %             uold_x_nu_new = (interp_local*[uold_x_nu(i);uold_x_nu(1)])';
        end
        ref_level_ele_new=[ref_level_ele, ref_level_ele] + 1;
    else
        if i<length(x)
            ele_new = [x(i),x(i+1)];
            uold_nu_new_ele = [uold_nu(i),uold_nu(i+1)];
            uold_x_nu_new_ele = [uold_x_nu(i),uold_x_nu(i+1)];
        elseif i==length(x)
            %             ele_new=  [x(i),Lx];
        end
        ref_level_ele_new = ref_level_ele;
    end
    
    if i>1
        ele_new=ele_new(2:end);
        uold_nu_new_ele=uold_nu_new_ele(2:end);
        uold_x_nu_new_ele=uold_x_nu_new_ele(2:end);
    end
    x_new = [x_new,ele_new];
    uold_nu_new = [uold_nu_new, uold_nu_new_ele];
    uold_x_nu_new = [uold_x_nu_new, uold_x_nu_new_ele];
    ref_level_new= [ref_level_new, ref_level_ele_new];
    
    
end

if size(x_new)~=size(ref_level_new)
    disp('disaster, sizes of x and ref_level not same')
end

end

function BD2 = fn_2S(dist_x_pl,uold_nu)
BD2 = (1./( ((dist_x_pl+ circshift(dist_x_pl,1))./dist_x_pl).^2 .*(-dist_x_pl) - (-dist_x_pl - circshift(dist_x_pl,1))   ))...
    .* (((dist_x_pl+ circshift(dist_x_pl,1))./dist_x_pl).^2 .*circshift(uold_nu,1) - circshift(uold_nu,2) - ...
    (((dist_x_pl+ circshift(dist_x_pl,1))./dist_x_pl).^2 -1).*uold_nu);
end

function f = fn_central_nonu(dist_x_pl, dist_x_min, u)
f = (1./dist_x_pl - 1./(dist_x_pl + dist_x_min)) .* circshift(u,-1)...
    +(1./dist_x_min - 1./(dist_x_pl)) .* u...
    +(1./(dist_x_pl+dist_x_min) - 1./(dist_x_min)) .*circshift(u,1);

% f= 1./(dist_x_pl+dist_x_min) .* (circshift(u,-1)+circshift(u,1))...
%     + (dist_x_pl -dist_x_min)./(dist_x_pl+dist_x_min)...
%     .*( ((circshift(u,-1) -u )./(dist_x_pl)) -  ((circshift(u,0) -circshift(u,1) )./(dist_x_min)));
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


