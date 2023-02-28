function var_out = fn_plot_lin_advect(u_init,v_init, ylim_plot1, ylim_plot2, xlim_plot, EI_plot_y_lim,R_plot_lims,pow,scheme_name_str, init_conds_str, interpolant_str)


c={};
scheme_arr = {scheme_name_str};
exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'pwl'};



for i = 1 : length(scheme_arr)
    c{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_',init_conds_str,'_',interpolant_str,'_lin_advect.mat']));
    %     c{i}=struct2cell('LeapFrog_cell_arr_file_sin_IC_uniform.mat');
    
end
% colour_arr =['y';'c';'g';"mo";"ro";"bo";"ko"];

colour_arr =['y';'c';'g';'m';'r';'b';'k'];
n_entries = length(c);
which_scheme = 1;

xlim_p = .1;
ylim_p = [0, 5];
n_schemes= length(c);
n_plots = 4;
f_size= 14;
N_subplots=4;
%% plots

i_int=1;
for i_scheme = 1:n_schemes
    %
   
        figgy= figure();
        
        l_refs = length(c{i_scheme}{1});
 for i_ref = 1: l_refs
        
        time_arr  =  c{i_scheme}{1}{i_ref}(1,:);
        bound_arr =  c{i_scheme}{1}{i_ref}(2,:);
        error_arr =  c{i_scheme}{1}{i_ref}(3,:);
        EI_arr    =  c{i_scheme}{1}{i_ref}(4,:);
        
        
        
        % evolution of the error
        subplot(1,5,1)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,error_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([0,xlim_p])
%         ylim([1e-11,1e-4])
        hold on;
        
        % Evolution of the indicator
        subplot(1,5,2)
        semilogy(time_arr,bound_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
        %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
        %         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        ylabel(['$\mathcal{E}\left(t^n\right)$'],'Interpreter','latex','FontSize',f_size)
        xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
        grid on;
        xlim([0,xlim_p])
%         ylim([1e-11,1e-4])
        hold on;
        
        
        if i_ref<l_refs
            time_arr_2  = c{i_scheme}{1}{i_ref+1}(1,:);
            bound_arr_2 = c{i_scheme}{1}{i_ref+1}(2,:);
            error_arr_2 = c{i_scheme}{1}{i_ref+1}(3,:);
            dofs_arr_2 = c{i_scheme}{1}{i_ref+1}(4,:);

            
            time_arr_1  = c{i_scheme}{1}{i_ref}(1,:);
            bound_arr_1 = c{i_scheme}{1}{i_ref}(2,:);
            error_arr_1 = c{i_scheme}{1}{i_ref}(3,:);
            dofs_arr_1 = c{i_scheme}{1}{i_ref}(4,:);

            %             length(time_arr_2(1:2:end))
            %             length(time_arr_1)
            
            if length(error_arr_2)>length(error_arr_1)
                EOC_error = log((error_arr_2(1:(2^(i_ref))^(1):end))./(error_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_bound = log((bound_arr_2(1:(2^(i_ref))^1:end))./(bound_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            else
                EOC_error = log((error_arr_2)./(error_arr_1))./log((dofs_arr_1)./(dofs_arr_2));
                EOC_bound = log((bound_arr_2)./(bound_arr_1))./log((dofs_arr_1)./(dofs_arr_2));
            end
            
            subplot(1,5,3)
            if length(time_arr_1)>length(EOC_error)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            %             ylabel(['$EOC\left( \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel(['$t^n$'],'Interpreter','latex','FontSize',f_size)
            grid on;
            %             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            xlim([0,xlim_p])
            ylim([0,2])
            hold on;
            
            subplot(1,5,4)
            if length(time_arr_1)>length(EOC_bound)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\mathcal{E}\left(t^n\right)\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([0,xlim_p])
%             ylim([0,6])
            hold on;
            %                 ylim([0,3])
            
            
            
            
        end
        
        
        
        
        subplot(1,5,5)
        plot(time_arr,EI_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
        
        ylabel(['$EI\left(t^n\right)$'],'Interpreter','latex','FontSize',f_size)
        xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
        xlim([0,xlim_p])
        ylim([0,10])
        %         ylim([0,10])
        
        grid on;
        hold on;
    end
        hold off;
        set(gcf, 'Position',  [100, 100, 1050, 200])
%          print(gcf,['/Users/gs1511/Desktop/GSialounas_repo_temporal_rec/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_harmonic_u',num2str(round(10*u_init)),'_v',num2str(round(10*v_init)),'_paperrec_poly_tristan.png'],'-dpng','-r600')
         print(gcf,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/figures/fig_',scheme_arr{i_scheme},'plots_1x5_step_IC_P1_shw'],'-dpng','-r600')

%             saveas(figgy,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_',scheme_arr{i_scheme},'_','plots_1x5_lin_advect'],'png');%char(Scheme_string)
%             saveas(figgy,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_uniform'],'png');%char(Scheme_string)
        
   
    
    %     sg%title([scheme_arr{i_scheme},' bound'],'Color','red')
    %     sgt.FontSize = 20;
    
end

hold off;

var_out =1;end