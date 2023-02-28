function var_out = fn_plot_lin_advect_L1_thesis(u_init,v_init, ylim_plot1, ylim_plot2, xlim_plot, EI_plot_y_lim,R_plot_lims,pow,scheme_name_str, init_conds_str, interpolant_str)

c={};
scheme_arr = {scheme_name_str};

for i = 1 : length(scheme_arr)
    c{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_disc_adv_for_thesis_',init_conds_str,'_',interpolant_str,'_lin_advect_comparison_adapt_OFF.mat']));
    %     c{i}=struct2cell('LeapFrog_cell_arr_file_sin_IC_uniform.mat');
    
end
% colour_arr =['y';'c';'g';"mo";"ro";"bo";"ko"];

colour_arr =['y';'c';'g';'m';'r';'b';'k'];
n_entries = length(c);

xlim_p = .5;
ylim_p = [-3, 3];
n_schemes= length(c);
n_plots = 4;
f_size= 18;
N_subplots=5;
%% plots


figgy= figure();

l_refs = length(c{1}{1});
for i_ref = 1: l_refs
    
    
    time_arr  =  c{1}{1}{i_ref}(1,:);
    bound_arr =  c{1}{1}{i_ref}(2,:);
    error_arr =  c{1}{1}{i_ref}(3,:);
    EI_arr    =  c{1}{1}{i_ref}(4,:);
    
    
    
    % evolution of the error
    subplot(1,N_subplots,1)
    semilogy(time_arr,error_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
    xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
    ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)

%     ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;\left|\left|\left|\cdot\right|\right|\right|\right)}$'],'Interpreter','latex','FontSize',f_size)
    grid on;
%     ylim([1e-11, 1e-5])
    xlim([0, xlim_p])
    %             ylim(ylim_p)
    hold on;
    
    % Evolution of the indicator
    subplot(1,N_subplots,2)
    semilogy(time_arr,bound_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
%     ylabel(['$\eta\left(t^n\right)$'],'Interpreter','latex','FontSize',f_size)
    ylabel(['$\left(\omega\mathcal{E}^2\right)^{1/2}$'],'Interpreter','latex','FontSize',f_size)

    xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
    grid on;
    xlim([0, xlim_p])
%     ylim([1e-11, 1e-5])
    hold on;
    
    if i_ref<l_refs
        time_arr_2  = c{1}{1}{i_ref+1}(1,:);
        bound_arr_2 = c{1}{1}{i_ref+1}(2,:);
        error_arr_2 = c{1}{1}{i_ref+1}(3,:);
        %
        %                 EOC_R1_arr_2 = c{i_scheme}{1}{i_ref+1}(6,:);
        %                 EOC_R2_arr_2 = c{i_scheme}{1}{i_ref+1}(7,:);
        %
        %                 EOC_R1_arr_1 = c{i_scheme}{1}{i_ref}(6,:);
        %                 EOC_R2_arr_1 = c{i_scheme}{1}{i_ref}(7,:);
        
        time_arr_1  = c{1}{1}{i_ref}(1,:);
        bound_arr_1 = c{1}{1}{i_ref}(2,:);
        error_arr_1 = c{1}{1}{i_ref}(3,:);
        
        if length(error_arr_2)>length(error_arr_1)
            EOC_error = log((error_arr_2(1:(2^(i_ref))^(1):end))./(error_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            EOC_bound = log((bound_arr_2(1:(2^(i_ref))^1:end))./(bound_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            
            %                     EOC_R1 = log((EOC_R1_arr_2(1:(2^(i_ref))^(1):end))./(EOC_R1_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            %                     EOC_R2 = log((EOC_R2_arr_2(1:(2^(i_ref))^(1):end))./(EOC_R2_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            
        else
            EOC_error = log((error_arr_2)./(error_arr_1))/log(0.5);
            EOC_bound = log((bound_arr_2)./(bound_arr_1))/log(0.5);
        end
        
        subplot(1,N_subplots,3)
        if length(time_arr_1)>length(EOC_error)
            plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
        else
            plot(time_arr_1,EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
        end
        
        %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        %             ylabel(['$EOC\left( \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
        ylabel(['$EOC\left(\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}\right)$'],'Interpreter','latex','FontSize',f_size)

%         ylabel(['$EOC\left(\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;\left|\left|\left|\cdot\right|\right|\right|\right)}\right)$'],'Interpreter','latex','FontSize',f_size)
        xlabel(['$t^n$'],'Interpreter','latex','FontSize',f_size)
        grid on;
        %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        xlim([0, xlim_p])
        ylim([0,4])
        hold on;
        
        subplot(1,N_subplots,4)
        if length(time_arr_1)>length(EOC_bound)
            plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
        else
            plot(time_arr_1,EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
        end
        ylabel(['$EOC\left(\eta\left(t^n\right)\right)$'],'Interpreter','latex','FontSize',f_size)
        xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
        %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        grid on;
        xlim([0, xlim_p])
        ylim([0,4])
        hold on;

        
    end
    
    
    
    
    subplot(1,N_subplots,5)
    %                         plot(time_arr,EI_arr+(l_refs-i_ref)*.05,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
    
    plot(time_arr(1:end),EI_arr(1:end),colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
    
    ylabel(['$EI\left(t^n\right)$'],'Interpreter','latex','FontSize',f_size)
    xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
    xlim([0, xlim_p])
%     ylim([0,EI_plot_y_lim])
    ylim([0,10])
    
    grid on;
    hold on;
    
 
    
end
hold off;
set(gcf, 'Position',  [100, 100, 1050, 200])
print(gcf,['/Users/gs1511/Desktop/GSialounas_PhD_Thesis_repo/figures/fig_linadvect_WENO_plots_1x5_L1_bound_',init_conds_str,'_.png'],'-dpng','-r600')

% saveas(figgy,['/Users/gs1511/Desktop/GSialounas_PhD_Thesis_repo/figures/fig_leapfrog_plots_1x5_u_wave_weno_rec'],'png');%char(Scheme_string)
%             saveas(figgy,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_uniform'],'png');%char(Scheme_string)


hold off;
var_out = 1;
end