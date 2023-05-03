clc;
clear;
close all;

global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();
% In this file we will run linear advection using SSP3 for temporal
% discretization and WENO5 for the following conditions:
% Initial conditions: smooth (sinusoidal), hat, step (discontinuous)
% Boudnary conditions: periodic:
% scheme_arr = {'SSP3WENO3_postshock'};
% scheme_arr = {'LxF_postshock'};
% scheme_str = 'SSP3WENO3_preshock';
% scheme_str = 'LxF_preshock';%{'LxF_preshock'};
% scheme_arr = {'LxW_postshock'};
% scheme_str = 'LxF_postshock';
scheme_str = 'FTBS';
% scheme_str = 'SSP3WENO3_postshock';

% error L1
% error L2
init_conds  = "stpIC";

% init_conds  = "pwlin";
% rec_string = 'weno';
% rec_string = 'hermite';
% rec_string= 'linear';
% rec_string = 'second';
for i=1:1
    
    if init_conds == "sinIC"
        N_pts_burger= 100;
        ex = @(x,t) burger_sol(x,t,N_pts_burger);%fn_burger_sol_lin_3(x,t);
        t_0  = 1e-100;
        init_cond = 'sinIC';
    elseif init_conds == "pwlin"
        init_cond = 'pwlin';
        ex = @ (x,t) fn_burger_ex_lin(x,t);
        t_0 =  0.5;
    elseif init_conds == "stpIC"
        init_cond = 'stpIC';
        ul = 0;
        ur = 1;
        ex = @(x,t) fn_exact_step_burger(x,t,ul,ur);
%         ex = @(x,t) fn_ex_transport_stp(x,t);

        t_0 = 0;
    end
    var_out = fn_ftbs_burger_fresh(scheme_str,init_cond,ex,t_0);
    
    fn_plot_FTBS_burger_L2L1_bound_comparison_offset(scheme_str,init_cond,t_0,ul,ur)
    
    
    
    pow_plot = 1;
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

