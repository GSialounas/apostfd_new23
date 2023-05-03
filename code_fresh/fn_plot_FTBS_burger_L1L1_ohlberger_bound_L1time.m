function var_out = fn_plot_FTBS_burger_L1L1_ohlberger_bound_L1time(scheme_str, init_conds,t_0,ul,ur)
c={};

exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'qdr'};

scheme_arr  ={scheme_str};

for i = 1 : length(scheme_arr)
    c{i}=struct2cell(load([scheme_str,'_cell_arr_file_',init_conds,'_burger_try_interp_ul_',num2str(ul),'_ur_',num2str(ur),'_ohlberger_bound.mat']));
end


colour_arr =['y';'c';'g';'m';'r';'b';'k'];
n_entries = length(c);
which_scheme = 1;

n_schemes= length(c);
n_plots = 4;
f_size= 14;
y_lim_error_L1 = [ 1e-4,1e-1];
y_lim_error_L2 = [1e-4,5e-1];
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
        bound_arr_l1 =  c{i_scheme}{1}{i_ref}(2,:);
        error_arr_l1  =  c{i_scheme}{1}{i_ref}(3,:); % Linfinity L1 error
        EI_arr =  c{i_scheme}{1}{i_ref}(4,:);
        ohlberger_arr = c{i_scheme}{1}{i_ref}(5,:);
        EI_arr_ohlberger =  c{i_scheme}{1}{i_ref}(6,:);
        error_arr_l1l1 =  c{i_scheme}{1}{i_ref}(7,:);
        error_arr_l2 =  c{i_scheme}{1}{i_ref}(8,:);
        bound_arr_l2 = c{i_scheme}{1}{i_ref}(9,:);
        min_bound = min(bound_arr_l2, ohlberger_arr);
        EI_arr_min  =min_bound./error_arr_l1l1;
        
        
        xlim_p = min(time_arr(end),2.2);
        % evolution of the error

        subplot(3,5,1)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,error_arr_l2,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim([1e-2 1e+1])
        hold on;
        
        
        subplot(3,5,6)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,error_arr_l1l1,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\left|\left|e\right|\right|_{L^{1}\left(0,t^n;L^1\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L1)
        hold on;
        
                
        subplot(3,5,11)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,error_arr_l1l1,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\left|\left|e\right|\right|_{L^{1}\left(0,t^n;L^1\left(\Omega\right)\right)}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L1)
        hold on;
        
        
        subplot(3,5,2)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,bound_arr_l2,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\left(\omega\mathcal{E}^2\left(t,L^2\left(\Omega\right)\right)\right)^{1/2}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim([1e-2 1e+1])
        hold on;
        
        subplot(3,5,7)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,ohlberger_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\omega\mathcal{E}\left(t,L^1\left(\Omega\right)\right)_{ohl}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L1)
        hold on;
%         
         
        subplot(3,5,12)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,min(ohlberger_arr,bound_arr_l2),[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$\omega\mathcal{E}\left(t,L^1\left(\Omega\right)\right)_{ohl}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
        ylim(y_lim_error_L1)
        hold on;
%         
     

        
        
        if i_ref<l_refs
            time_arr_2  = c{i_scheme}{1}{i_ref+1}(1,:);
            bound_arr_l1_2 = c{i_scheme}{1}{i_ref+1}(2,:);
            error_arr_l1_2 = c{i_scheme}{1}{i_ref+1}(3,:);
            ohlberger_arr_l1_2 = c{i_scheme}{1}{i_ref+1}(5,:);
            error_arr_l1l1_2 = c{i_scheme}{1}{i_ref+1}(7,:);
            error_arr_l2_2 = c{i_scheme}{1}{i_ref+1}(8,:);
            bound_arr_l2_2 = c{i_scheme}{1}{i_ref+1}(9,:);

            min_arr_2 = min(bound_arr_l2_2, ohlberger_arr_l1_2);
            
            time_arr_1  = c{i_scheme}{1}{i_ref}(1,:);
            bound_arr_l1_1 = c{i_scheme}{1}{i_ref}(2,:);
            error_arr_l1_1 = c{i_scheme}{1}{i_ref}(3,:);
            ohlberger_arr_l1_1 = c{i_scheme}{1}{i_ref}(5,:);
            error_arr_l1l1_1 = c{i_scheme}{1}{i_ref}(7,:);
            error_arr_l2_1 = c{i_scheme}{1}{i_ref}(8,:);
            bound_arr_l2_1 = c{i_scheme}{1}{i_ref}(9,:);

            min_arr_1 = min(bound_arr_l2_1, ohlberger_arr_l1_1);

            
            %             length(time_arr_2(1:2:end))
            %             length(time_arr_1)
            
            if length(bound_arr_l1_2)>length(bound_arr_l1_1)
 
                EOC_bound_l1 = log((bound_arr_l1_2(1:(2^(i_ref))^1:end))./(bound_arr_l1_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_error_l1 = log((error_arr_l1_2(1:(2^(i_ref))^1:end))./(error_arr_l1_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_ohlberger_l1 = log((ohlberger_arr_l1_2(1:(2^(i_ref))^1:end))./(ohlberger_arr_l1_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_error_l1l1 = log((error_arr_l1l1_2(1:(2^(i_ref))^1:end))./(error_arr_l1l1_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_bound_l2 = log((bound_arr_l2_2(1:(2^(i_ref))^1:end))./(bound_arr_l2_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_error_l2 = log((error_arr_l2_2(1:(2^(i_ref))^1:end))./(error_arr_l2_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_min_arr = log((min_arr_2(1:(2^(i_ref))^1:end))./(min_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            else
                              
                EOC_bound_l1 = log((bound_arr_l1_2)./(bound_arr_l1_1))/log(0.5);
                EOC_ohlberger_l1 = log((ohlberger_arr_l1_2)./(ohlberger_arr_l1_1))/log(0.5);
                EOC_error_l1 = log((error_arr_l1_2)./(error_arr_l1_1))/log(0.5);
                EOC_error_l1l1 = log((error_arr_l1l1_2)./(error_arr_l1l1_1))/log(0.5);
                EOC_bound_l2 = log((bound_arr_l2_2)./(bound_arr_l2_1))/log(0.5);
                EOC_error_l2 = log((error_arr_l2_2)./(error_arr_l2_1))/log(0.5);
                EOC_min_arr = log((min_arr_2)./(min_arr_1))/log(0.5);


            end
         
            
             subplot(3,5,3)
            if length(time_arr_1)>length(EOC_bound_l1)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error_l2,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_error_l2,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\left|\left|e\right|\right|_{L^{\infty}\left(0,t^n;L^2\left(\Omega\right)\right)}\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
            
            subplot(3,5,8)
            if length(time_arr_1)>length(EOC_bound_l1)
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
            
            subplot(3,5,13)
            if length(time_arr_1)>length(EOC_bound_l1)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error_l1l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_error_l1l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(min\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
            subplot(3,5,4)
            if length(time_arr_1)>length(EOC_bound_l1)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_bound_l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_bound_l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\left(\omega\mathcal{E}^2\left(t,L^2\left(\Omega\right)\right)\right)^{1/2}\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
            subplot(3,5,9)
            if length(time_arr_1)>length(EOC_bound_l1)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_ohlberger_l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_ohlberger_l1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(\omega\mathcal{E}\left(t,L^1\left(\Omega\right)\right)\right)_{ohl}$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
             subplot(3,5,14)
            if length(time_arr_1)>length(EOC_bound_l1)
                plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_min_arr,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_min_arr,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end
            
            %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(min\right)$'],'Interpreter','latex','FontSize',f_size)
            xlabel('$t^n$','Interpreter','latex','FontSize',f_size)
            %             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([time_arr(1),xlim_p])
            ylim(y_lim_EOC)
            hold on;
            
           
        end

         subplot(3,5,5)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,EI_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$EI\left(t^n\right)$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
%         ylim(y_lim_error_L1)
        hold on;
        
        subplot(3,5,10)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,EI_arr_ohlberger,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$EI\left(t^n\right)_{ohl}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
%         ylim(y_lim_error_L1)
        hold on;
        
         subplot(3,5,15)
        %         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        semilogy(time_arr,EI_arr_min,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t^n$','Interpreter','Latex','FontSize',f_size)
        ylabel(['$EI\left(t^n\right)_{min}$'],'Interpreter','latex','FontSize',f_size)
        
        grid on;
        xlim([time_arr(1),xlim_p])
%         ylim(y_lim_error_L1)
        hold on;
       
    end
    
    hold off;
    set(gcf, 'Position',  [100, 100, 1050, 550])
%     saveas(figgy,['/Users/gs1511/Desktop/GSialounas_PhD_Thesis_repo/figures/fig_SSP3WENO3_preshock',init_conds,'_plots_1x5_burgers'],'png');%char(Scheme_string)
    if ul>ur
        condition = 'shock';
    else
        condition = 'rarefaction';
    end
    
    if unique(init_conds == 'hatIC')==1
        condition = 'hatIC';
    elseif unique(init_conds == 'cmbIC')==1
        condition = 'cmbIC';
    else
        
    end
    print(gcf,['/Users/gs1511/Desktop/apostfd_new23/figures/fig_',scheme_str,'_bound_comparison_burgers_offset_',num2str(10*t_0),'_',condition,'_simpleinterp_ohlberger_bound_L1_time'],'-dpng','-r600')
    saveas(figgy,['/Users/gs1511/Desktop/apostfd_new23/figures/fig_',scheme_str,'_bound_comparison_burgers_offset_',num2str(10*t_0),'_',condition,'_simpleinterp_low_res_ohlberger_bound_L1_time'],'png');%char(Scheme_string)

    
    hold off;
    var_out =1;
end