clc;
clear;
close all;

record_video =0;
rep_video=1;
mesh_adapt=0;
if mesh_adapt == 1
    adaptivity = 'adapt_ON';
else
    adaptivity = 'adapt_OFF';
end

global nq xq wq
nq = 2;
gauss();

% T=4*pi;
showplots = 1;
save_results_to_cellarr=1;
L_domain= 32*pi;
time_stepping = 'RK3';
rec ='rec_1';
spat_disc = 'WENO3';
dt_dx_coupling = 'fixed';
scheme_string= ['SHW_dambreak',time_stepping,'_',spat_disc,'_',rec,'_',dt_dx_coupling,'_gs_',adaptivity];
scheme_arr = {scheme_string};

fps =1;
g=9.81;

h_l = .2;
h_r = .1;
x_dam = 30;
[c_m, h_m] = fn_shw_dambreak_cm_hm(g, h_l, h_r);
fn_exact_soln = @(x,t) fn_shw_dambreak_exact(t,x, x_dam, L_domain, h_l, h_r,g, c_m,h_m);
fluxfun='shallow_water'; % select flux function
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
    case 'shallow_water'
        % w = [h; hv]
        flux = @(w) [w(2,:); (w(2,:).^2./w(1,:))+.5*g*w(1,:).^2];
        eig_Jf_1 = @(w) w(2,:)./w(1,:) + sqrt(g*w(1,:));
        eig_Jf_2 = @(w) w(2,:)./w(1,:) - sqrt(g*w(1,:));
        dflux = @(w) [0 , 1; -(w(2,:).^2./(w(1,:).^2)) + g*w(1,:), 2*w(2)./w(1)];
end


sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        %         S = @(w) 0.1*w.^2;
        
        %          S = @(w) .1*w.^2;
        %         S =@(x,t) fn_exact_forcing(x,t,d,a,omega,k,g);
    case 'dont'
        S = @(w) zeros(size(w));
end

% scheme_arr = {'SHW_RK1_LXF_rec1_pow_two'};

x = linspace(0,L_domain,251);
dist_x_pl = x(2:end)- x(1:end-1);
% dt = .05 *(x(2)-x(1));
% the bcs are non-periodic so we consider the last grid point, x=L, as
% well.
% x = x(1:end-1);


f_size = 14;
N_subplots = 3;
p=0; q=L_domain;


i_exponent = 1;
i_interpolant = 1;
i_scheme =1;
cell_cell_arr_shw = {};

flux_fn  = @(dist_x_pl,dt,Ut) ones(size(Ut));

maxit_arr = [1:4]; % max refinement iterations
axis_p = [0 q -3e-3 3e-3];
dt_arr = L_domain*(2.^(-(maxit_arr+7))); % mesh size
dx_arr = sqrt(10*dt_arr);

for m = 1:1
    t= 0;
    
    it =1;
    
    dx = L_domain*(2^(-(maxit_arr(m)+8))); % mesh size
    dx_coarsest = L_domain*(2^(-(maxit_arr(1)+3)));
    dx_finest   = L_domain*(2^(-(maxit_arr(end)+8)));
    dt_coarsest = .05* dx_coarsest;
    
    
    % Choice of dt;
    if dt_dx_coupling == "fixed"
        dt=.1*dx;
    elseif dt_dx_coupling=="pow_one"
        dt=.5*dx;
    elseif dt_dx_coupling=="pow_two"
        dt=.1*dx_finest^2;
    else
    end
    
    
    if m ==1
        T= 3000*dt;
    end
    %     x=0:dx:L_domain;
    ratio_cTf = 1;
    dx_fine= dx;
    part_fine = 1;
    x = create_grid(L_domain, ratio_cTf, dx_fine, part_fine)   ;
    dist_x_pl = x(2:end)-x(1:end-1);
    dist_x_min = circshift(dist_x_pl,1);
    %     x = x(1:end-1);
    %     h= 5+cos(2*pi*x);
    %     v=sin(cos(2*pi*x))./h;
    %%initial conditions
    %periodic
    %     h = d+a*sin(k*x);
    %     v = a*omega*(cosh(k*h)/sinh(k*d)).*sin(omega*t-k*x);
    %     dam-break
    h = fn_hinit_dambreak(x);
    v = zeros(size(h));
    h_old = [h;h.*v];
    
    
    L2L2R = 0; %L2L2 accumulation of ||R||^2
    L2L2R_h = 0; %L2L2 accumulation of ||R_h||^2
    L2L2R_hv = 0; %L2L2 accumulation of ||R_hv||^2
    L2L2R_arr = zeros(1, length(x)-1);
    
    
    % arrays used for EOC for bound and error
    bound_arr = []; %zeros(1,ceil(T/dt));
    error_arr = []; % zeros(1,ceil(T/dt));
    time_arr  = []; % zeros(1,ceil(T/dt));
    dofs_arr  = []; % zeros(1,ceil(T/dt));
    EI_index  = [];
    
    time_arr(1) = 0;
    
    %
    % form initial mesh
    tree_list = [];
    root_ids = [];
    
    tol = 1e-6;
    max_ref= 4;
    for i = 1:length(x)-1
        if i==1
            tree_list = [ Tree_vector(element_vector(x(i), x(i+1),i,tol),i,tol)];
            %                 tree_list(i) = tree_list(i).assign_vals([sin(x(i)), sin(x(i+1))], exp(-100*(x(i)-.5).^2));
            root_ids(i) = i;
        else
            tree_list(end+1) = Tree_vector(element_vector(x(i), x(i+1),i,tol),i,tol);
            %                 tree_list(i) = tree_list(i).assign_vals([sin(x(i)), sin(x(i+1))], exp(-100*(x(i)-.5).^2));
            root_ids(i) = i;
        end
    end
    
    if mesh_adapt==1
        M = mesh_vector(tree_list,root_ids,max_ref);
        M.updateMesh();
        %         x = x(1:end-1);
        M.asgnVals(h_old, zeros(size(h_old)));
        
        M.refineGlobal();
        M.refineGlobal();
        M.refineGlobal();
        M.refineGlobal();
%         M.refineGlobal();
%         M.refineGlobal();



        [x_ref,h_old_ref] = M.getNewVals();
        
        % amend the grid and uld
        %         x_ref = [x_ref, Lx];
        h_old  = h_old_ref;
        
        dist_x_pl_ref = x_ref(2:end) -x_ref(1:end-1);
        dist_x_min_ref = circshift(dist_x_pl_ref, 1);
        dist_x_pl = dist_x_pl_ref;
        dist_x_min = dist_x_min_ref;
        x=x_ref;
        %         x = x_ref(1:end-1);
    else
        %         x=x(1:end-1);
    end
    if mesh_adapt == 1
        h= fn_hinit_dambreak(x);
        v = zeros(size(h));
        M.asgnVals([h;v], zeros(size([h;v])));
    end
    if mesh_adapt==1
    dt = 2*2*.1*(x(2)-x(1)); % four ref levels
    T = 2*2500*dt;
    
%     dt = .1*(x(2)-x(1)); % five ref levels
%     T = 12*2500*dt;
    else
        dt= .1 *(x(2)-x(1));
        T= 3000*dt;
    end
    while  t < (T-dt/2)
        [un1, un2, F1, F2] = dependent(h,v,g);
%         f_h_old = [flux_lxf_dambreak(dist_x_pl,dt,un1,F1);flux_lxf_dambreak(dist_x_pl,dt,un2,F2)];

        h_old = [un1; un2];
        hq_old = h_old;
%         f_h_old_weno = WENO3resAdv1d_FD_system_just_for_rec_dambreak([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2,x,t);
        f_h_old_weno = WENO3_FD_system_dambreak_nonu([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2,x,t);

        %         h_new = h_old - dt*f_h_old;
%         h = h_new(1,:);
%         v= h_new(2,:)./h_new(1,:);
%         [un1, un2, F1, F2] = dependent(h,v,g);
%         f_h_new = [flux_lxf_dambreak(dist_x_pl,dt,un1,F1);flux_lxf_dambreak(dist_x_pl,dt,un2,F2)];

        if time_stepping== "RK3"
            hq_stage1 = hq_old - dt* f_h_old_weno;
            h_stage1 = hq_stage1(1,:);
            v_stage1 = hq_stage1(2,:)./h_stage1;
            [un1_stage1, un2_stage1, F1_stage1, F2_stage1] = dependent(h_stage1,v_stage1,g);
            
%             f_h_stage1_weno = WENO3resAdv1d_FD_system_just_for_rec_dambreak([un1_stage1;un2_stage1],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1_stage1,F2_stage1,x,t);
            f_h_stage1_weno = WENO3_FD_system_dambreak_nonu([un1_stage1;un2_stage1],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1_stage1,F2_stage1,x,t);

            % Stage 2
            hq_stage2 = .75* hq_old +.25*hq_stage1 -.25*dt*f_h_stage1_weno;
            h_stage2 = hq_stage2(1,:);
            v_stage2 = hq_stage2(2,:)./h_stage2;
            
            [un1_stage2, un2_stage2, F1_stage2, F2_stage2] = dependent(h_stage2,v_stage2,g);
%             f_h_stage2_weno = WENO3resAdv1d_FD_system_just_for_rec_dambreak([un1_stage2;un2_stage2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1_stage2,F2_stage2,x,t);
            f_h_stage2_weno = WENO3_FD_system_dambreak_nonu([un1_stage2;un2_stage2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1_stage2,F2_stage2,x,t);

            % Stage 3
            hq_stage3 = (1/3)* hq_old +(2/3)*hq_stage2 -(2/3)*dt*f_h_stage2_weno;
            h_stage3 = hq_stage3(1,:);
            v_stage3 = hq_stage3(2,:)./h_stage3;
            
            h_new = [hq_stage3(1,:);hq_stage3(2,:)];
            h = h_new(1,:);
            v = h_new(2,:)./h;
            
            [un1, un2, F1, F2] = dependent(h,v,g);
            
%             f_h_new_weno = WENO3resAdv1d_FD_system_just_for_rec_dambreak([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2,x,t);
            f_h_new_weno = WENO3_FD_system_dambreak_nonu([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2,x,t);

        end
        
        if mesh_adapt==0
            L2L2R_arr= zeros(1,length(x)-1);
            L2L2R_arr_h= zeros(1,length(x)-1);
            L2L2R_arr_hv= zeros(1,length(x)-1);
            for iq = 1 : nq
                tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                if rec =="rec_1"
                    [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv, IUt_h, IUt_hv,RL2iq_h_arr,RL2iq_hv_arr,] = compute_Rs_vector_temp_1_spatiotemp_1_shw_dambreak(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'WENO3',flux_fn); %compute R at gauss point
                elseif rec =="rec_3"
                    % this gives second order until the shock appears but has a
                    % central flux
                    %                 [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv,IUt_h, IUt_hv] = compute_Rs_vector_temp_3_spatiotemp_3_shw(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'CS',flux_fn); %compute R at gauss point
                    [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, IUt_h, IUx_hv,IUt_hv, IUx_hv_flux ] = compute_Rs_vector_temp_3_spatiotemp_3_shw_weno3(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'CS',flux_fn,flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2); %compute R at gauss point
                    
                else
                end
                L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                L2L2R_arr_h = L2L2R_h + wq(iq)*dt*(RL2iq_h); %quadrature formula
                L2L2R_arr_hv = L2L2R_hv + wq(iq)*dt*(RL2iq_hv); %quadrature formula
                L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
                
                if (L2L2R<0)
                    disp('L2L2R<0')
                end
            end
        else
            L2L2R_arr= zeros(1,length(x)-1);
            L2L2R_arr_h= zeros(1,length(x)-1);
            L2L2R_arr_hv= zeros(1,length(x)-1);
            for iq = 1 : nq
                tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                if rec =="rec_1"
                    [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv, IUt_h, IUt_hv,RL2iq_h_arr,RL2iq_hv_arr,] = compute_Rs_vector_temp_1_spatiotemp_1_shw_dambreak(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old,f_h_new,'CS',flux_fn); %compute R at gauss point
                elseif rec =="rec_3"
                    % this gives second order until the shock appears but has a
                    % central flux
                    %                 [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv,IUt_h, IUt_hv] = compute_Rs_vector_temp_3_spatiotemp_3_shw(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'CS',flux_fn); %compute R at gauss point
                    [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, IUt_h, IUx_hv,IUt_hv, IUx_hv_flux ] = compute_Rs_vector_temp_3_spatiotemp_3_shw_weno3(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'CS',flux_fn,flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2); %compute R at gauss point
                    
                else
                end
                %                 L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                L2L2R_arr_h = L2L2R_arr_h + wq(iq)*dt*(RL2iq_h); %quadrature formula
                L2L2R_arr_hv = L2L2R_arr_hv + wq(iq)*dt*(RL2iq_hv); %quadrature formula
                L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
                
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
            h = h_old(1,:);
            v= h_old(2,:)./h_old(1,:);
            M.asgnVals([h;v], RL2iq_arr);
            M.mark_ref_coar();
            M.ref_coarsen_local();
            [x_ref, h_old_ref] = M.getNewVals();
            
            % amend the grid and uld
            %             x_ref = [x_ref, Lx];
            h  = h_old_ref(1,:);
            v  = h_old_ref(2,:);
            
            dist_x_pl_ref = x_ref(2:end) -x_ref(1:end-1);
            dist_x_min_ref = circshift(dist_x_pl_ref, 1,2);
            dist_x_pl = dist_x_pl_ref;
            dist_x_min = dist_x_min_ref;
            x = x_ref;
            
            [un1, un2, F1, F2] = dependent(h,v,g);
            h_old = [un1; un2];
            f_h_old = [1./dist_x_pl; 1./dist_x_pl] .* [F1(2:end)- F1(1:end-1); F2(2:end) - F2(1:end-1)];
            [f_h_new, h_new] = fn_lxw_dambreak_modified(un1,un2./un1,dist_x_pl, dt,g);
            
            % take the time-step with the refined uold
            
            %                 dt = .1 *min(dist_x_pl);
%             [un1, un2, F1, F2] = dependent(h,v,g);
%             h_old = [un1;un2];
%             hq_old = h_old;
            %         f_h_old_weno = WENO3resAdv1d_FD_system([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2,x,t);
%             f_h_old = [flux_lxf_dambreak(dist_x_pl,dt,un1,F1);flux_lxf_dambreak(dist_x_pl,dt,un2,F2)];
            %         f_h_old  = WENO5resAdv1d_fdm_gs(uold,flux,dflux,S,dx);%resWENO5(uold_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu); %WENO5resAdv1d_fdm_gs(uold_nu,flux,dflux,S,dx);%
            
            %         f_h_old = flux_richtmayer_h(dist_x_pl,dt,un1, un2,F1, F2,g);
            
            %         un = [un1;un2] - dt* f_h_old;
            %         forcing_term = S(x,t);
            % stage 1:
            %             hq_stage1 = hq_old - dt* f_h_old_weno;
%             hq_stage1 = hq_old - dt* f_h_old;
            
%             h_new = [hq_stage1(1,:);hq_stage1(2,:)];
            h = h_new(1,:);
            v = h_new(2,:)./h;
            [un1, un2, F1, F2] = dependent(h,v,g);
            %             f_h_new_weno = WENO5resAdv1d_fdm_gs_system([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2,x,t);
%             f_h_new=[flux_lxf_dambreak(dist_x_pl,dt,un1,F1);flux_lxf_dambreak(dist_x_pl,dt,un2,F2)];
            % bound for plotting purposes
            
            L2L2R_arr= zeros(1,length(x)-1);
            L2L2R_arr_h= zeros(1,length(x)-1);
            L2L2R_arr_hv= zeros(1,length(x)-1);
            for iq = 1 : nq
                tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                if rec =="rec_1"
                    [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv, IUt_h, IUt_hv,RL2iq_h_arr,RL2iq_hv_arr,] = compute_Rs_vector_temp_1_spatiotemp_1_shw_dambreak(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old,f_h_new,'CS',flux_fn); %compute R at gauss point
                elseif rec =="rec_3"
                    % this gives second order until the shock appears but has a
                    % central flux
                    %                 [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv,IUt_h, IUt_hv] = compute_Rs_vector_temp_3_spatiotemp_3_shw(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'CS',flux_fn); %compute R at gauss point
                    [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, IUt_h, IUx_hv,IUt_hv, IUx_hv_flux ] = compute_Rs_vector_temp_3_spatiotemp_3_shw_weno3(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'CS',flux_fn,flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2); %compute R at gauss point
                    
                else
                end
                L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                L2L2R_arr_h = L2L2R_arr_h + wq(iq)*dt*(RL2iq_h); %quadrature formula
                L2L2R_arr_hv = L2L2R_arr_hv + wq(iq)*dt*(RL2iq_hv); %quadrature formula
                L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
                
                if (L2L2R<0)
                    disp('L2L2R<0')
                end
            end
        end
        
        if it ==0 && mesh_adapt==0
            c_0_coeff_arr_0 = c_0_coeff_arr_old;
            x_0 = x;
            dist_x_pl_0 = dist_x_pl;
            dist_x_min_0 = dist_x_min;
        end
        
        %
        %         if it ==0
        %             % t=0;
        %             bound_arr(1) = sqrt(exp(it*dt)*(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
        %             error_new = sqrt(space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,u,0*dt,ex,c_0_coeff_arr_0));
        %             error_arr(1) =error_new;% max(error_new, error_old);
        %             error_old = max(error_arr);
        %             EI_index(1) = bound_arr(1)./error_arr(1);
        %             error_new_arr(end+1) = error_new;
        %             dofs_arr(1) = length(x_0);
        %
        %             % t=dt
        %             bound_arr(end+1) = sqrt(exp(1*dt)*(L2L2R + space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
        %             error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,1*dt,ex,c_0_coeff_arr_new));
        %             error_arr(end+1) =error_new;% max(error_new, error_old);
        %             error_old = max(error_arr);
        %             EI_index(end+1) = bound_arr(end)./error_arr(end);
        %             error_new_arr(end+1) = error_new;
        %             dofs_arr(end+1) = length(x);
        %
        %
        %         else
        %             % t>dt
        %             bound_arr(end+1) = sqrt(exp((it+1)*dt)*(L2L2R + space_int_vector_gs(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),0,ex,c_0_coeff_arr_0)));
        %             error_new = sqrt(space_int_vector_gs(x,dist_x_pl,dist_x_min,u,(it+1)*dt,ex,c_0_coeff_arr_new));
        %             error_arr(end+1) = error_new;%max(error_new, error_old);
        %             error_old = max(error_arr);
        %             EI_index(end+1) = bound_arr(end)./error_arr(end);
        %             error_new_arr(end+1) = error_new;
        %             dofs_arr(end+1) = length(x);
        %
        %         end
        
        time_arr(it) = (it-1)*dt;
        bound_arr_eoc(it) = sqrt(L2L2R);
        bound_arr_h_eoc(it) = sqrt(L2L2R_h);
        bound_arr_hv_eoc(it) = sqrt(L2L2R_hv);
        
        % error_arr_eoc(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,it*dt,ex));
        error_arr_eoc(it) = 1;%sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,it*dt,ex,c_0_coeff_arr_new));
        dofs_arr(it) = length(x);
        EI_index(it) = bound_arr_eoc(it)./error_arr_eoc(it);
        it = it + 1;
        t = t + dt  ;
        
        if (mod(it-1,fps)==0 &&showplots==1)
            %             S_forcing = S(x,t);
            
            [ex_soln_h, ex_soln_hv] = fn_exact_soln(x,t);
%             ex_soln_h = ex_soln(1,:);
%             ex_soln_hv = ex_soln(2,:);
            subplot(2,2,1)
            grid on;
            
            %             plot(x,h,'b',x, ex_soln_h,'ro')
            plot(x,h,'bo',x,ex_soln_h,'r')
            ylabel('$h$','Interpreter','latex','FontSize',f_size)
                xlabel('$x$','Interpreter','latex','FontSize',f_size)
            grid on;
            axis([0 L_domain .1 .21])
            title(t)
            
            subplot(2,2,2)
            plot(x(1:end-1),RL2iq_arr,'r')
            ylabel(['$\mathcal{E}\left(t^n; \left[x_j, x_{j+1}\right]\right)$'],'Interpreter','latex','FontSize',f_size)

            grid on;
            %             plot(x,S_forcing(1,:),'b')
            ylim([0, 1e-5])
            %             axis([])
            %             title('IU_h')
            title('R_h')
            
            
            
         
            subplot(2,2,3)
            
            plot(x,h.*v,'bo', x, ex_soln_hv,'r')
            ylabel('$hv$','Interpreter','latex','FontSize',f_size)
                xlabel('$x$','Interpreter','latex','FontSize',f_size)
            %             plot(x,v,'b')
            axis([0 L_domain 0 .1])
            
            grid on;
            
            %             plot(x,h.*v,'b', x, ex_soln_hv,'ro')
            
            %             axis([0 L_domain 0 .2])
            %             title('hv')
            subplot(2,2,4)
            
            plot(x(1:end-1),dist_x_pl,'b.')
             ylabel(['$Grid\, Spacing$'],'Interpreter','latex','FontSize',f_size)
                xlabel('$x$','Interpreter','latex','FontSize',f_size)
            grid on;
            %             axis(axis_p)
            %             title('IU_hv')
            title('R_{hv}')
            
            pause(0.01)
            set(gcf, 'Position',  [100, 100, 800, 400])
            
            if (record_video)
                F_var(rep_video) = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                
            end
            rep_video = rep_video+1;
            %             if ((m==length(maxit_arr))&&(reps==20|| reps==1700))
            %                 saveas(figgy_parasite,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_',scheme_arr{i_scheme},'_',num2str(reps),'_dtdx_',dt_dx_coupling,'_',num2str(m)],'png')
            %             end
        end
        
        h = h_new(1,:);
        v= h_new(2,:)./h(1,:);
        
    end
    
    if (record_video)
            
            video = VideoWriter(['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/videos/WENO3_shw_dambreak_ic','_',num2str(maxit_arr(m)),'_solution_with_',adaptivity,'.mp4'],'MPEG-4');
            video.FrameRate = fps;
            open(video);
            writeVideo(video,F_var);
            close(video);
        end
    if (save_results_to_cellarr==1)
        cell_cell_arr_shw{m}=[time_arr;bound_arr_eoc;error_arr_eoc;EI_index;dofs_arr;bound_arr_h_eoc;bound_arr_hv_eoc];
    end
end
save([scheme_arr{i_scheme},'_cell_arr_file_shw_',adaptivity,'.mat'],'cell_cell_arr_shw')




%
% function f_h = flux_richtmayer(dist_x_pl,dt,u,F)
% u_pl  = .5 * (circshift(u,-1) + u) - (dt./(2*dist_x_pl)) .* (circshift(F,-1) - F);
% u_min = .5 * (circshift(u,1) + u) - (dt./(2*dist_x_pl)) .* (F-  circshift(F,1));
% f_h = (1./dist_x_pl) .* .5 .* (u_pl.^2  - u_min.^2);
% end
%
%
function f_out = flux_richtmayer_h(dist_x_pl,dt,h, hv,F_h, F_hv,g)
h_pl  = .5 * (circshift(h,-1) + h) - (dt./(2*dist_x_pl)) .* (circshift(F_h,-1) - F_h);
h_min = .5 * (circshift(h,1) + h) - (dt./(2*dist_x_pl)) .* (F_h-  circshift(F_h,1));

hv_pl  = .5 * (circshift(hv,-1) + hv) - (dt./(2*dist_x_pl)) .* (circshift(F_hv,-1) - F_hv);
hv_min = .5 * (circshift(hv,1) + hv) - (dt./(2*dist_x_pl)) .* (F_hv-  circshift(F_hv,1));

f_h_pl = hv_pl;
f_h_min = hv_min;

f_hv_pl = (hv_pl.^2)./(h_pl) +  .5* g * h_pl.^2;
f_hv_min = (hv_min.^2)./(h_min) +  .5* g * h_min.^2;


f_h  = (1./dist_x_pl) .* (f_h_pl-f_h_min);
f_hv = (1./dist_x_pl) .* (f_hv_pl-f_hv_min);

f_out= [f_h; f_hv];
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
