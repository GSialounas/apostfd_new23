clc;
clear;
close all;


record_video =1;
rep_video=1;
mesh_adapt=1;
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
time_stepping = 'RK1';
rec ='rec_1';
spat_disc = 'LxW';
dt_dx_coupling = 'fixed';
scheme_string= ['SHW_',time_stepping,'_',spat_disc,'_',rec,'_',dt_dx_coupling,'_gs_',adaptivity];
scheme_arr = {scheme_string};

fps =10;
g=9.81;
d= 2;
lambda = 32*pi;
k = 2*pi/lambda;
omega = sqrt(g*k*tanh(k*d));
a=.1;
h_l = .2;
h_r = .1;
x_dam = 30;
[c_m, h_m] = fn_shw_dambreak_cm_hm(g, h_l, h_r);
fn_exact_soln = @(x,t) fn_shw_dambreak_exact(t,x, x_dam, L_domain, h_l, h_r,g, c_m,h_m);



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
    
    ratio_cTf = 1;
    dx_fine= dx;
    part_fine = 1;
    x = create_grid(L_domain, ratio_cTf, dx_fine, part_fine)   ;
    dist_x_pl = x(2:end)-x(1:end-1);
    dist_x_min = circshift(dist_x_pl,1);
    
    x_mid_eval = x(2:end)- x(1:end-1);
    dist_x_pl_mid_eval = x_mid_eval(2:end)-x_mid_eval(1:end-1);
    dist_x_min_mid_eval = circshift(dist_x_pl_mid_eval,1);
    
    h = fn_hinit_dambreak(x_mid_eval);
    v = zeros(size(h));
    h_old = [h;h.*v];
    
end