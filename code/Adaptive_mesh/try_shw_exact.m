clc;
clear;
close all;

L_domain = 32*pi;
x= linspace(0,L_domain,501);


h = fn_hinit_dambreak(x);
dt = .1*(x(2)-x(1));
x_dam = 30;
h_l = .2;
h_r = .1;
v = zeros(size(x));
g= 9.81;


coeff_1 = 1;
coeff_2 = 0;
coeff_3 = -9*g*h_r;
coeff_4 = 16*g*h_r*sqrt(g*h_l);
coeff_5 = -8*g^2*h_l*h_r - g^2*h_r^2;
coeff_6 = 0;
coeff_7 = g^3*h_r^3;
p = [coeff_1, coeff_2, coeff_3, coeff_4, coeff_5, coeff_6, coeff_7];
r= roots(p);
r_criterion= r.^2/g;

h_m =-5;
for i =1:length(r_criterion)
    if real(r_criterion(i))>=h_r && real(r_criterion(i))<h_l
        h_m = r_criterion(i);
    end
end

T= 25;
t= dt;
h=zeros(size(x));
v = zeros(size(x));
figure
while t<T
    c_m =sqrt(g*h_m);
    x_A = x_dam - t*sqrt(g*h_l);
    x_B = x_dam+ t*(2*sqrt(g*h_l) -3*c_m);
    x_C = x_dam + t*(2*c_m^2 * (sqrt(g*h_l)-c_m))/(c_m^2 -g*h_r);
    
    h(x<=x_A) = h_l;
    h(x>x_A & x<= x_B) = (4/(9*g)) *(sqrt(g*h_l) - (x(x>x_A & x<= x_B)-x_dam)/(2*t)).^2;
    h(x>x_B & x<= x_C) = c_m^2/g;
    h(x>x_C) = h_r;
    
    v(x<=x_A) = 0;
    v(x>x_A & x<= x_B) = (2/3) *((x(x>x_A & x<= x_B)-x_dam)/t + sqrt(g*h_l)) ;
    v(x>x_B & x<= x_C) = 2 * (sqrt(g*h_l) - c_m);
    v(x>x_C) = 0;
    
    subplot(2,1,1)
    plot(x,h,'b')
    subplot(2,1,2)
    plot(x,h.*v,'b')
    pause(0.01)
    t= t+dt;
end

