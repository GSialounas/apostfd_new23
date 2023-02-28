clc;
clear;
close all;

showplots = 1;
L_domain= 32*pi;
g = 9.81; 

x = linspace(0,L_domain, 501);
dx = x(2) - x(1);
dt = .1*dx;

T = 10000*dt;
t = 0;

h= fn_hinit_dambreak(x);
v = zeros(size(x));

h = [0 h 0];
v = [ 0 v 0];

[h,v] = fn_boundary(h,v);
un1_new = zeros(size(h));
un2_new = zeros(size(h));

un1_pl = zeros(size(h));
un2_pl = zeros(size(h));

un1_min = zeros(size(h));
un2_min = zeros(size(h));
it =0;
fps = 10;
[un1, un2, F1, F2] = dependent(h,v,g); 
figure
while t<(T-dt/2)
    
    for j = 2:length(h)-1
        un1_pl(j) = .5 * (un1(j+1)+un1(j)) - .5 * (dt/dx) * (F1(j+1) - F1(j));
        un2_pl(j) = .5 * (un2(j+1)+un2(j)) - .5 * (dt/dx) * (F2(j+1) - F2(j));
        
        un1_min(j) = .5 * (un1(j)+un1(j-1)) - .5 * (dt/dx) * (F1(j) - F1(j-1));
        un2_min(j) = .5 * (un2(j)+un2(j-1)) - .5 * (dt/dx) * (F2(j) - F2(j-1));
        
    end
    % Calculate plus values
    h = un1_pl;
    v=  un2_pl./un1_pl;
    [h,v] = fn_boundary(h,v);
    [~,~,F1_pl, F2_pl] = dependent(h,v,g);
    
    % Calculate minus values
    h = un1_min;
    v=  un2_min./un1_min;
    [h,v] = fn_boundary(h,v);
    [~,~,F1_min, F2_min] = dependent(h,v,g);
    
    for j = 2:length(h) -1
        un1_new(j) = un1(j) - (dt/dx) * (F1_pl(j) - F1_min(j));
        un2_new(j) = un2(j) - (dt/dx) * (F2_pl(j) - F2_min(j));
    end
    
    h = un1_new;
    v= un2_new./un1_new;
    [h,v] = fn_boundary(h,v);
    [un1, un2, F1, F2] = dependent(h,v,g);
    t= t+dt;
    it = it+1;
    subplot(2,1,1)
    plot(x,h(2:end-1),'b',x, fn_hinit_dambreak(x),'r')
    title(t)
    subplot(2,1,2)
    plot(x,v(2:end-1),'r')
    pause(0.01)
    
    
    
end