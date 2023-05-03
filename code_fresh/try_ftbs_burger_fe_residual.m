clc;
clear;
close all;

fps = 10;
showplot = 1;
global mesh uh fh nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();
scheme_str = 'FTBS';
init_cond = 'stpIC';
ul = 0;
ur = 1;

if init_cond == 'stpIC'
    ex_step = @(x,t) fn_exact_step_burger(x,t,ul,ur);
    ex_step_x = @(x,t) zeros(size(x));
else
    ex_smooth  = @(x,t) sin(10*pi*x);
    ex_smooth_x = @(x,t) 10*pi*cos(10*pi*x);
end


cell_arr = {};

for m = 1:4
    disp(m)
    h = 2^(-m-6);
    dt= .1*h;
    if m ==1
        T = 800*dt;
    end
    x = -2:h:2;
    dist_x_pl =  x(2:end) -x(1:end-1);
    it =0;
    t = 0;
    
    %% Create FE mesh
    n_vertices = length(x);
    n_elem = n_vertices -1;
    for i = 1:length(x)
        if i <length(x)
            elem_vertices(i,:) = [i,i+1];
        end
        vertex_coordinates(i,:) = x(i);
    end
    
    mesh.n_elem = n_elem;
    mesh.n_vertices= n_vertices;
    mesh.elem_vertices = elem_vertices;
    mesh.vertex_coordinates =  vertex_coordinates;
    
    elem_neighbours=  zeros(n_elem,2);
    elem_boundaries = zeros(n_elem,2);
    elem_boundaries(1,1)= 1;
    elem_boundaries(end,end)=1;
    mesh.elem_boundaries = elem_boundaries;
    for i = 1:length(elem_neighbours)
        if i ==1
            elem_neighbours(i,2) = 2;
        elseif i == length(elem_neighbours)
            elem_neighbours(i,1) = n_elem-1;
        else
            elem_neighbours(i,:) = [i-1,i+1];
        end
    end
    mesh.elem_neighbours= elem_neighbours;
    uh = zeros(mesh.n_vertices,1);
    fh = zeros(mesh.n_vertices,1);
    
    ndof(m) = mesh.n_vertices;
    
    
    
    % ic


    
    L1L1R = 0;
    h_arr(m) = h;
    x=x(1:end-1);
    mask  = ones(size(x));
    mask(((x-x(1))<=4*h) | ((x(end)-x)<=4*h))=0;
    
    N= length(x);
    dx = x(2)-x(1);
    S=  h * spdiags([1/6*ones(N,1), 2/3*ones(N,1), 1/6*ones(N,1)], -1:1, N,N);
    %      B=B;
    S(1,:)=0;
    S(end,:)=0;
    S(1,1) = 1;
    S(end,end)=1;
    
    
    uold = ex_step(x,t);
    u = uold;
    bound_arr = [];
    time_arr= [];
    error_0 = 0;
    
    for iq = 1:nq
        xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
        diff_x = xiq-x;
        IU = fn_WENO3_rec(xiq, x,dist_x_pl, uold);        
        error_0 = error_0 + sum(wq(iq).* dist_x_pl .* abs(IU - ex_step(xiq,0)));
    end
    bound_arr(1) = error_0;
    time_arr(1) = 0;
    
    figure
    
    while length(time_arr)<2^(m-1)*800+1
%         u = uold - dt/h*.5*(uold.^2 - circshift(uold,1).^2);
        u = uold - dt/h*(uold.^1 - circshift(uold,1).^1);
        
        L1L1R_arr = zeros(size(x));
        for iq = 1:nq
            tq(iq) = 0.5*dt*xq(iq) + dt/2 +t;
%             [R1L1_iq, R1L1_iq_arr] = compute_temp_1_space_weno(x,dist_x_pl,uold,u,tq(iq),t,dt);
            [R_i,R1L1_iq,R1L1_iq_arr]  = assemble_and_solve_1d_fe_residual(S,x,uold,u,tq(iq),t,dt,ul,ur);
          
            L1L1R = L1L1R + dt*wq(iq)*R1L1_iq;
            L1L1R_arr = L1L1R_arr + dt*wq(iq)*R1L1_iq_arr;
        end
        
        if (showplot==1 && mod(it,fps)==0)
            
            plot(x,u,'b')
            ylim([-.1 1.1])
%             xlim([-1 1])
            pause(0.01)
        end
        t = t+dt;
        it = it+1;
        uold = u;
        bound_arr(end+1) = L1L1R + error_0;
        time_arr(end+1) = it*dt;
        
    end
    
     
    cell_arr{m} = [time_arr;bound_arr];
end
    save(['FTBS_cell_arr_file_',init_cond,'_burger_try_interp_fe_residual_ul_',num2str(ul),'_ur_',num2str(ur),'.mat'],'cell_arr')
    t_0=0;
    
    fn_plot_FTBS_burger_L1L1_bound_fe_residual(scheme_str,init_cond,t_0,ul,ur)

colour_arr =['y';'c';'g';'m';'r';'b';'k'];
    figure
    hold on;
    l_refs=4;
    for i =1:4
        time_arr = cell_arr{i}(1,:);
        bound_arr =  cell_arr{i}(2,:);
        semilogy(time_arr,bound_arr,[colour_arr(end-(l_refs-i))])
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