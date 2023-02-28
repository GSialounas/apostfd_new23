function [c_m, h_m] = fn_shw_dambreak_cm_hm(g, h_l, h_r)
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
c_m =sqrt(g*h_m);
end