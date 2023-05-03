clc;
clear;
close all;

global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();

init_cond = 'stpIC';
ul = 1;
ur = 0;
ex_step = @(x,t) fn_exact_step_burger(x,t,ul,ur);
ex_step_x = @(x,t) zeros(size(x));
ex_smooth  = @(x,t) sin(10*pi*x);
ex_smooth_x = @(x,t) 10*pi*cos(10*pi*x);

ex_hat = @(x,t) fn_hat(x);

evalt= 0 ;
N = 6;
error_arr = zeros(1,N);
error_arr_derivative = zeros(1,N);
h_arr = zeros(1,N);
err_weno_old = -1;
err_weno_derivative_old = -1;
EOC_arr = zeros(size(error_arr));
EOC_arr_derivative = zeros(size(error_arr));
for m = 1:length(error_arr)
    h = 2^(-m-3);
    x = -2:h:2;
    
    L2  = 0 ;
    h_arr(m) = h;
    x=x(1:end-1);
    mask  = ones(size(x));
    mask(((x-x(1))<=4*h) | ((x(end)-x)<=4*h))=0;
    
    for iq = 1:nq
        xiq = x +.5*xq(iq)*h + h/2;
        [IU, IU_x] = fn_WENO3_rec(xiq,x,h,ex_step(x,0));
        error_arr(m) = error_arr(m) + sum(wq(iq)*h.*mask.*abs(IU - ex_step(xiq,evalt)));
        error_arr_derivative(m) = error_arr_derivative(m) + sum(wq(iq)*h*mask.*abs(IU_x - ex_step_x(xiq,evalt)));
    end
    
    figure
    subplot(1,2,1)
    plot(xiq,IU,'b',xiq,ex_step(xiq,evalt),'r')
    title('IU vs ex')
    subplot(1,2,2)
    plot(xiq,IU_x,'b',xiq,ex_step_x(xiq,evalt),'r')
    title('IU_x vs ex_x')
    
    if err_weno_old <0
        EOC_err(m)=0;
        EOC_err_derivative(m) = 0;
    else
        EOC_arr(m) = log(err_weno_old/error_arr(m))/log(2);
        EOC_arr_derivative(m) = log(err_weno_derivative_old/error_arr_derivative(m))/log(2);
    end
    err_weno_old = error_arr(m);
    err_weno_derivative_old = error_arr_derivative(m);
end

F_size = 14;
figure
subplot(2,2,1)
loglog(1./h_arr,error_arr,'b--o')
xlabel('$h^{-1}$','Interpreter','latex')
ylabel('$\left|\left| u-\widehat{U} \right|\right|_{\mathrm{L}^1\left(\Omega\right)}$','Interpreter','latex')

subplot(2,2,2)
loglog(1./h_arr,error_arr_derivative,'b--o')
xlabel('$h^{-1}$','Interpreter','latex')
ylabel('$\left|\left| u_x-\widehat{U}_x \right|\right|_{\mathrm{L}^1\left(\Omega\right)}$','Interpreter','latex')


subplot(2,2,3)
semilogx(1./h_arr(2:end),EOC_arr(2:end),'k:+','Linewidth',1)
% hold on;
% loglog(ndof_arr,err_arr_hermite,'k:^','Linewidth',1)
% hold on;
% loglog(ndof_arr,ndof_arr.^-2,'k--','Linewidth',1)
xlabel('$h^{-1}$','Interpreter','latex','FontSize',F_size)
ylabel('$EOC\left(||u-\widehat{U}||_{\mathrm{L}^2\left(\Omega\right)}\right)$','Interpreter','latex','FontSize',F_size)
% legend('$\left|\left| \nabla e\right|\right|_{\mathrm{L}^2}$','$\sum \left(\eta_j^2\right)^{1/2}$','Interpreter','latex','Location','SouthWest','Fontsize',14)
ylim([0 6])
grid on;
hold off;
% 
subplot(2,2,4);
semilogx(1./(h_arr(2:end)),EOC_arr_derivative(2:end),'k:+','Linewidth',1)
ylabel('$EOC\left(||u_x-\widehat{U}_x||_{\mathrm{L}^2\left(\Omega\right)}\right)$','Interpreter','latex','FontSize',F_size)

xlabel('$h^{-1}$','Interpreter','latex','FontSize',F_size);
ylim([0 6])
grid on;



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