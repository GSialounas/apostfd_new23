clc;
clear;
close all;

x = 0:0.01:20;
y1 = 200*exp(-0.05*x).*sin(x);
y2 = 200*exp(-0.05*x).*cos(x);
% [AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
% set(AX,{'ycolor'},{'b';11*[0.05 0.05 0.05]})  % Left color red, right color blue...
f_size = 14;
%
yyaxis right
plot(x, y2,'LineWidth',1,'Color',16*[0.05 0.05 0.05])

set(gca,'ycolor',11*[0.05 0.05 0.05]) 
ylabel('$\ln(1/h_j)$','Interpreter','latex','FontSize',f_size)

yyaxis left
h=plot(x,y1,'LineWidth',1)
%             plot(x,h,x,ex_soln_h)
% xlim([0 32*pi])
% ylim([0.05,.25])
%
% %             plot(x,h-ex_soln_h,'r')
% hold on;
ylabel('$\eta$','Interpreter','latex','FontSize',f_size)
xlabel('$x_j$','Interpreter','latex','FontSize',f_size)
grid on;

