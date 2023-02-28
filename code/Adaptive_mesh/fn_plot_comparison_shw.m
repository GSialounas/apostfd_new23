function var_out = fn_plot_comparison_shw(u_init,v_init, ylim_plot1, ylim_plot2, xlim_plot, EI_plot_y_lim,R_plot_lims,pow,scheme_name_str, init_conds_str, interpolant_str)


c_ON={};
c_OFF={};

scheme_arr = {scheme_name_str};
exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'pwl'};



for i = 1 : length(scheme_arr)
    c_ON{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_',init_conds_str,'_',interpolant_str,'_lin_advect_comparison_adapt_ON.mat']));
    c_OFF{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_',init_conds_str,'_',interpolant_str,'_lin_advect_comparison_adapt_OFF.mat']));

    %     c{i}=struct2cell('LeapFrog_cell_arr_file_sin_IC_uniform.mat');
    
end
colour_arr_OFF =['y';'c';'g';"mo";"ro";"b--o";"k"];

colour_arr =['y';'c';'g';'m';'r';'b';"r"];
n_entries = length(c_ON);
which_scheme = 1;

xlim_p = .5;
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
        
        
        
        
        % evolution of the error
%         subplot(1,N_subplots,1)
%         %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%         semilogy(time_arr_ON,error_arr_ON,[colour_arr(end-(l_refs-i_ref))],time_arr_OFF,error_arr_OFF,[colour_arr_OFF(end-(l_refs-i_ref))],'LineWidth',2)
%         xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
% %         ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
%         ylabel(['$\left|\left|e\left(t^n\right)\right|\right|_{L^2\left(\Omega\right)}$'],'Interpreter','latex','FontSize',f_size)
%  
%         grid on;
%         xlim([0,xlim_p])
%         ylim([1e-2,2e-1])
%         hold on;
%         
%         subplot(1,N_subplots,2)
%         %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%         semilogy(time_arr_ON,bound_arr_ON,[colour_arr(end-(l_refs-i_ref))],time_arr_OFF,bound_arr_OFF,[colour_arr_OFF(end-(l_refs-i_ref))],'LineWidth',2)
%         xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
% %         ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
%         ylabel(['$\left|\left|e\left(t^n\right)\right|\right|_{L^2\left(\Omega\right)}$'],'Interpreter','latex','FontSize',f_size)
%  
%         grid on;
%         xlim([0,xlim_p])
% %         ylim([1e-2,2e-1])
%         hold on;
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
%         ylim([1e-2,2e-1])
        hold on;
        
         subplot(1,N_subplots,2)
        plot(cumsum(dofs_arr_ON),error_arr_ON,colour_arr(end-(l_refs-i_ref)),cumsum(dofs_arr_OFF),error_arr_OFF,[colour_arr_OFF(end-(l_refs-i_ref))],'LineWidth',2)
        %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
        %         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        ylabel(['$\left|\left|e\left(t^n\right)\right|\right|_{L^2\left(\Omega\right)}$'],'Interpreter','latex','FontSize',f_size)
        xlabel(['$\sum_{n=0}^N N_{dofs}\left(t^n\right)$'],'Interpreter','latex','FontSize',f_size)
        grid on;
        xlim([0, max(cumsum(dofs_arr_ON))])
%         xlim([0,xlim_p])
%         ylim([0,1100])
        hold on;
        legend('adaptive','uniform','Location','SouthEast')
        
%          subplot(1,3,3)
%         %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%         semilogy(time_arr_ON,bound_arr_ON,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
%         xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
% %         ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
%         ylabel(['$\left|\left|e\left(t^n\right)\right|\right|_{L^2\left(\Omega\right)}$'],'Interpreter','latex','FontSize',f_size)
%  
%         grid on;
%         xlim([0,xlim_p])
% %         ylim([1e-2,2e-1])
%         hold on;
        
%         
%          % evolution of the error
%         subplot(2,2,3)
%         %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%         semilogy(time_arr_OFF,error_arr_OFF,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
%         xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
%         ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
%         
%         grid on;
%         xlim([0,xlim_p])
%         ylim([1e-2,8e-2])
%         hold on;
%         
%         % Evolution of the indicator
%         subplot(2,2,4)
%         plot(time_arr_OFF,dofs_arr_OFF,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
%         %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
%         %         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
%         ylabel(['$N\left(dofs\right)$'],'Interpreter','latex','FontSize',f_size)
%         xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
%         grid on;
%         xlim([0,xlim_p])
%         ylim([0,4500])
%         hold on;
    
    end
        hold off;
        set(gcf, 'Position',  [100, 100, 800, 200])
%          print(gcf,['/Users/gs1511/Desktop/GSialounas_repo_temporal_rec/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_harmonic_u',num2str(round(10*u_init)),'_v',num2str(round(10*v_init)),'_paperrec_poly_tristan.png'],'-dpng','-r600')
         print(gcf,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/figures/fig_',scheme_arr{i_scheme},'plots_1x5_step_IC_P1_comparison_adaptONOFF'],'-dpng','-r600')
        
         arr_wr_names = ["time_arr_ON", "error_arr_ON", "error_arr_OFF"];
         arr_wr = [time_arr_ON', error_arr_ON', error_arr_OFF'];
         fid  = fopen('../../paper/data/error_file_adapt.csv','w');
         fprintf(fid,'%s, %s, %s\n',arr_wr_names(1), arr_wr_names(2), arr_wr_names(3))
         
         for i =1:length(arr_wr)
             fprintf(fid, '%5.4f,%5.4f,%5.4f\n', arr_wr(i,:));
         end
         fclose(fid);
         
         arr_wr_names = ["time_arr_ON", "dofs_arr_ON", "dofs_arr_OFF"];
         arr_wr = [time_arr_ON', dofs_arr_ON', dofs_arr_OFF'];
         fid  = fopen('../../paper/data/dofs_file_adapt.csv','w');
         fprintf(fid,'%s, %s, %s\n',arr_wr_names(1), arr_wr_names(2), arr_wr_names(3))
         
         for i =1:length(arr_wr)
             fprintf(fid, '%5.4f,%5.4f,%5.4f\n', arr_wr(i,:));
         end
         fclose(fid);
%          file = [time_arr_ON',error_arr_ON'];
%          xlswrite('/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/figures/My_file.xlsx',file); 

%             saveas(figgy,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_',scheme_arr{i_scheme},'_','plots_1x5_lin_advect'],'png');%char(Scheme_string)
%             saveas(figgy,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_uniform'],'png');%char(Scheme_string)
        
   
    
    %     sg%title([scheme_arr{i_scheme},' bound'],'Color','red')
    %     sgt.FontSize = 20;
    
end

hold off;

var_out =1;end