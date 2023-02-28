function [h, hv] = fn_shw_dambreak_exact(t,x, x_dam, L_domain, h_l, h_r,g, c_m,h_m)



h=zeros(size(x));
v = zeros(size(x));


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

hv= h.*v;
end