function var_out = fn_plot_comparison_shw_WENO_dambreak_thesis_corrections(u_init,v_init, ylim_plot1, ylim_plot2, xlim_plot, EI_plot_y_lim,R_plot_lims,pow,scheme_name_str, init_conds_str, interpolant_str)


c_ON={};
c_OFF={};

scheme_arr = {scheme_name_str};
exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'pwl'};



for i = 1 : length(scheme_arr)
    c_ON{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_shw_adapt_ON.mat']));
    c_OFF{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_shw_adapt_OFF.mat']));
    

    
end
colour_arr_OFF =['y';'c';'g';"mo";"ro";"b--o";"k"];
colour_arr =['y';'c';'g';'m';'r';'b';"r"];
n_entries = length(c_ON);
which_scheme = 1;

xlim_p = inf;
ylim_p = [0, 5];
n_schemes= length(c_ON);
n_plots = 4;
f_size= 18;
N_subplots=2;
%% plots

i_int=1;
for i_scheme = 1:n_schemes
    %
   
        figgy= figure();
        
        l_refs = length(c_ON{i_scheme}{1});
 for i_ref = l_refs: l_refs
        
        time_arr_ON  =  c_ON{i_scheme}{1}{i_ref}(1,:);
        bound_arr_ON =  c_ON{i_scheme}{1}{i_ref}(2,:);
        error_arr_ON =  c_ON{i_scheme}{1}{i_ref}(3,:);
%         EI_arr    =  c_ON{i_scheme}{1}{i_ref}(4,:);
        dofs_arr_ON  =  c_ON{i_scheme}{1}{i_ref}(5,:);
        
        time_arr_OFF  =  c_OFF{i_scheme}{1}{i_ref}(1,:);
        bound_arr_OFF =  c_OFF{i_scheme}{1}{i_ref}(2,:);
        error_arr_OFF =  c_OFF{i_scheme}{1}{i_ref}(3,:);
%         EI_arr    =  c_OFF{i_scheme}{1}{i_ref}(4,:);
        dofs_arr_OFF  =  c_OFF{i_scheme}{1}{i_ref}(5,:);
        
        
    
        % get the L_inf(L2 norm of errors)
        max_ON=error_arr_ON(1);
        max_OFF= error_arr_OFF(1);
     
        error_arr_ON_max = zeros(size(error_arr_ON));
        error_arr_OFF_max = zeros(size(error_arr_OFF));
        for i = 1: length(time_arr_ON)
           error_arr_ON_max(i) = max(error_arr_ON(1:i));
           error_arr_OFF_max(i) = max(error_arr_OFF(1:i));
          
        end
        error_arr_ON = error_arr_ON_max;
        error_arr_OFF = error_arr_OFF_max;
       
%         
%         
        % Evolution of the indicator
        subplot(1,N_subplots,1)
        plot(time_arr_ON,dofs_arr_ON,colour_arr(end-(l_refs-i_ref)),time_arr_OFF,dofs_arr_OFF,[colour_arr_OFF(end-(l_refs-i_ref))],'LineWidth',2)
        %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
        %         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        ylabel(['$N_{dofs}\left(t^n\right)$'],'Interpreter','latex','FontSize',f_size)
        xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
        grid on;
        xlim([0,xlim_p])
        ylim([3e+2,5e+3])
        hold on;
        
         subplot(1,N_subplots,2)
        plot(cumsum(dofs_arr_ON),error_arr_ON_max,colour_arr(end-(l_refs-i_ref)),cumsum(dofs_arr_OFF),error_arr_OFF_max,[colour_arr_OFF(end-(l_refs-i_ref))],'LineWidth',2)
        %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
        %         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        ylabel(['$\left|\left|e\left(t^n\right)\right|\right|_{L^\infty\left(0;t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
        xlabel(['$\sum_{n=0}^N N_{dofs}\left(t^n\right)$'],'Interpreter','latex','FontSize',f_size)
        grid on;
        xlim([0, max(cumsum(dofs_arr_ON))])
%         xlim([0,xlim_p])
%         ylim([0,1100])
        hold on;
        legend('adaptive','uniform-equivalent','Location','SouthEast')
        on=cumsum(dofs_arr_ON);
        off = cumsum(dofs_arr_OFF);
           time_arr_ON  =  time_arr_ON(1:40:end);
        bound_arr_ON = bound_arr_ON(1:40:end);
        error_arr_ON =  error_arr_ON(1:40:end);
%         EI_arr    =  c_ON{i_scheme}{1}{i_ref}(4,:);
        dofs_arr_ON  =   dofs_arr_ON (1:40:end)
        
        time_arr_OFF  =  time_arr_OFF(1:40:end);
        bound_arr_OFF = bound_arr_OFF(1:40:end);
        error_arr_OFF =  error_arr_OFF(1:40:end);
%         EI_arr    =  c_OFF{i_scheme}{1}{i_ref}(4,:);
        dofs_arr_OFF  =  dofs_arr_OFF(1:40:end);
        
        
     
         
    
 end
        % take only 1 in every 100 points
        hold off;
        set(gcf, 'Position',  [100, 100, 800, 300])
%          print(gcf,['/Users/gs1511/Desktop/GSialounas_repo_temporal_rec/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_harmonic_u',num2str(round(10*u_init)),'_v',num2str(round(10*v_init)),'_paperrec_poly_tristan.png'],'-dpng','-r600')
         print(gcf,['/Users/gs1511/Desktop/GSialounas_PhD_Thesis_repo/figures/fig_',scheme_arr{i_scheme},'plots_1x5_shw_dambreak_comparison_adaptONOFF'],'-dpng','-r600')
      
    
end

hold off;

var_out =1;
end