clc;
clear;
close all;


x = linspace(-2,2,501);
dx = x(2)-x(1);
dt= .1*dx;
fps=10;
t= 0 ;
it=0;
figure;
while t<1.5
    if mod(it,fps)==0
        plot(x,fn_burger_ex_lin(x,t),'b')
   pause(0.01)
    end
    t=t+dt;
    it =it+1;
end