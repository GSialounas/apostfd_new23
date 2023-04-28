function var_out = fn_plot_FTBS_burger_L1L1_bound_comparison(scheme_str, init_conds,t_0)

c={};

exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'qdr'};

scheme_arr  ={scheme_str};

for i = 1 : length(scheme_arr)
    c{i}=struct2cell(load([scheme_str,'_cell_arr_file_',init_conds,'_burger_l1l1_bound_comparison.mat']));
end


colour_arr =['y';'c';'g';'m';'r';'b';'k'];
n_entries = length(c);
which_scheme = 1;

n_schemes= length(c);
n_plots = 4;
f_size= 10;
y_lim_error_L1 = [ 1e-4,1e+0];
y_lim_error_L2 = [1e-4,1e+0];
y_lim_EOC = [0,1.1];
ylim_EI = [0,2]
%% plots
i_int=1;
for i_scheme = 1:n_schemes
    figgy= figure();
    
    l_int = length(interpolant_arr); % we will only present dt~dx^2
    %     for  i_int = 1:l_int
    l_refs = length(c{i_scheme}{1});
    for i_ref = 1: l_refs
        
        time_arr  =  c{i_scheme}{1}{i_ref}(1,:);
        bound_arr =  c{i_scheme}{1}{i_ref}(2,:);
        error_arr =  c{i_scheme}{1}{i_ref}(3,:);
        bound_arr_l1 =  c{i_scheme}{1}{i_ref}(4,:);
        error_arr_l1 =  c{i_scheme}{1}{i_ref}(5,:);
        bound_arr_ohl =  c{i_scheme}{1}{i_ref}(6,:);
        error_arr_l1l1 =  c{i_scheme}{1}{i_ref}(7,:);
        
        
        
%         EI_arr    =  c{i_scheme}{1}{i_ref}(6,:);
        
        
        xlim_p = time_arr(end);
        % evolution of the error
        subplot(3,4,1)
        
        semilogy(time_arr,error_arr_l1,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^1\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L1)
        hold on;
        
        subplot(3,4,5)
        
        semilogy(time_arr,error_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L2)
        hold on;
        
        
        subplot(3,4,9)
        
        semilogy(time_arr,error_arr_l1l1,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\left|\left|e\right|\right|_{L^{1}\left(0,t^n;L^1\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L2)
        hold on;

        
        subplot(3,4,2)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,bound_arr_l1,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\omega\mathcal{E}\left(t,L^1\left(\Omega\right)\right)$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L1)
        hold on;
%         

        subplot(3,4,6)
        semilogy(time_arr,bound_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
        %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
        %         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        ylabel(['$\left(\omega\mathcal{E}^2\left(t,L^2\left(\Omega\right)\right)\right)^{1/2}$'],'Interpreter','latex','FontSize',f_size)
        xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L2)
        hold on;
        % Evolution of the indicator
        
        subplot(3,4,10)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,bound_arr_ohl,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\omega\mathcal{E}\left(t,L^1\left(\Omega\right)\right)_{ohl}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L1)
        hold on;

        
        
        if i_ref<l_refs
            time_arr_2  = c{i_scheme}{1}{i_ref+1}(1,:);
            bound_arr_2 = c{i_scheme}{1}{i_ref+1}(2,:);
            error_arr_2 = c{i_scheme}{1}{i_ref+1}(3,:);
            bound_arr_l1_2 = c{i_scheme}{1}{i_ref+1}(4,:);
            error_arr_l1_2 = c{i_scheme}{1}{i_ref+1}(5,:);
            bound_arr_ohl_2 = c{i_scheme}{1}{i_ref+1}(6,:);
            error_arr_l1l1_2 = c{i_scheme}{1}{i_ref+1}(7,:);
            
            time_arr_1  = c{i_scheme}{1}{i_ref}(1,:);
            bound_arr_1 = c{i_scheme}{1}{i_ref}(2,:);
            error_arr_1 = c{i_scheme}{1}{i_ref}(3,:);
            bound_arr_l1_1 = c{i_scheme}{1}{i_ref}(4,:);
            error_arr_l1_1 = c{i_scheme}{1}{i_ref}(5,:);
            bound_arr_ohl_1 = c{i_scheme}{1}{i_ref}(6,:);
            error_arr_l1l1_1 = c{i_scheme}{1}{i_ref}(7,:);
            %             length(time_arr_2(1:2:end))
            %             length(time_arr_1)
            
            if length(error_arr_2)>length(error_arr_1)
                EOC_error = log((error_arr_2(1:(2^(i_ref))^(1):end))./(error_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_bound = log((bound_arr_2(1:(2^(i_ref))^1:end))./(bound_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                
                EOC_error_l1 = log((error_arr_l1_2(1:(2^(i_ref))^(1):end))./(error_arr_l1_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_bound_l1 = log((bound_arr_l1_2(1:(2^(i_ref))^1:end))./(bound_arr_l1_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                
                EOC_error_l1l1 = log((error_arr_l1l1_2(1:(2^(i_ref))^(1):end))./(error_arr_l1l1_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_bound_ohl = log((bound_arr_ohl_2(1:(2^(i_ref))^1:end))./(bound_arr_ohl_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
 
            else
                EOC_error = log((error_arr_2)./(error_arr_1))/log(0.5);
                EOC_bound = log((bound_arr_2)./(bound_arr_1))/log(0.5);
                
                EOC_error_l1 = log((error_arr_l1_2)./(error_arr_l1_1))/log(0.5);
                EOC_bound_l1 = log((bound_arr_l1_2)./(bound_arr_l1_1))/log(0.5);
                
                EOC_error_l1l1 = log((error_arr_l1l1_2)./(error_arr_l1l1_1))/log(0.5);
                EOC_bound_ohl = log((bound_arr_ohl_2)./(bound_arr_ohl_1))/log(0.5);
            end
            
            
            subplot(3,4,3)
            if length(time_arr_1)>length(EOC_error)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error_l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_error_l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            %             ylabel(['$EOC\left( \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^1\left(\Omega\right)\right)}\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel(['$t^n$'],'Interpreter','latex','FontSize',f_size)
            grid on;
            %             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
            subplot(3,4,7)
            if length(time_arr_1)>length(EOC_bound)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
            subplot(3,4,11)
            if length(time_arr_1)>length(EOC_bound_ohl)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error_l1l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_error_l1l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\left|\left|e\right|\right|_{L^{1}\left(0,t^n;L^1\left(\Omega\right)\right)}\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
            
            subplot(3,4,4)
            if length(time_arr_1)>length(EOC_error)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_bound_l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_bound_l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            %             ylabel(['$EOC\left( \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\omega\mathcal{E}\left(t,L^1\left(\Omega\right)\right)\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel(['$t^n$'],'Interpreter','latex','FontSize',f_size)
            grid on;
            %             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
            subplot(3,4,8)
            if length(time_arr_1)>length(EOC_bound)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\left(\omega\mathcal{E}^2\left(t,L^2\left(\Omega\right)\right)\right)^{1/2}\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
            
            subplot(3,4,12)
            if length(time_arr_1)>length(EOC_error)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_bound_ohl,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_bound_ohl,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            %             ylabel(['$EOC\left( \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\omega\mathcal{E}\left(t,L^1\left(\Omega\right)\right)_{ohl}\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel(['$t^n$'],'Interpreter','latex','FontSize',f_size)
            grid on;
            %             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
           
        end
        
        
       
%         subplot(1,3,3)
%         semilogy(time_arr,bound_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
%         %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
%         %         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
%         ylabel(['$\left(\omega\mathcal{E}^2\left(t^n\right)\right)^{1/2}$'],'Interpreter','latex','FontSize',f_size)
%         xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
%         grid on;
%         xlim([0.98,1.02])
%         ylim([1e-1,1e+1])
%         hold on;
%         
%         
        
        
       
    end
    %     end
    
    %     sgtitle([scheme_arr{i_scheme},' bound'],'Color','red')
    %     sgt.FontSize = 20;
    hold off;
    set(gcf, 'Position',  [100, 100, 1050, 400])
%     saveas(figgy,['/Users/gs1511/Desktop/GSialounas_PhD_Thesis_repo/figures/fig_SSP3WENO3_preshock',init_conds,'_plots_1x5_burgers'],'png');%char(Scheme_string)
    print(gcf,['/Users/gs1511/Desktop/apostfd_new23/figures/fig_',scheme_str,'_bound_comparison_l1l1_burgers_offset_',num2str(10*t_0),'_',init_conds],'-dpng','-r600')

    
    hold off;
    var_out =1;
end