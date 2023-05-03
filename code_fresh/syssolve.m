clc;
clear;
close all;

A = [1 ,1 ; .25*(1+sqrt(5)), .25*(1-sqrt(5))];
x = [1;1];

m = A\x;

A1 = m(1);
A2 = m(2);


for n = 0:4
    rn = A1*(.25*(1+5^.5))^n + A2*(.25*(1-5^.5))^n;
    fprintf("n=%d and r_n = %d\n", [n,rn])
end

r_n = (sqrt(5)+3)/(2*sqrt(5)) *((1+sqrt(5))/4)^4 + (sqrt(5)-3)/(2*sqrt(5)) *((1-sqrt(5))/4)^4 ;

arr = 2:1002;

arr_sol  = arr.*((sqrt(5)+1)/4).^(arr);

plot(arr,cumsum(arr_sol))


arr = 1:1000;

m1 = arr./(2.^(arr+1));
m2 = (arr-1)./(2.^arr);

figure
plot(arr,cumsum(m1),'b', arr,cumsum(m2),'r')


m1 = (1+sqrt(5))/4;
m2 = (1-sqrt(5))/4;

sum1 = (m1/(m1-1))*m1^-2*(sqrt(5)+3)/(8*sqrt(5))*(4*2*m1/(3-sqrt(5)));
sum2 = (m2/(m2-1))*m2^-2*(sqrt(5)-3)/(8*sqrt(5))*(4*2*m2/(3+sqrt(5)));

sum_tot = sum1+sum2;

% sum_tot = fn_sum(1000);

factor_1 = (1/4) * (sqrt(5)+3)/(2*sqrt(5)) * ((1+sqrt(5))/4)^(-2)  ;
factor_2 = (1/4) * (sqrt(5)-3)/(2*sqrt(5)) * ((1-sqrt(5))/4)^(-2);



%% calculation

m1 = 4/(1+sqrt(5));
c1 = (sqrt(5)+3)/(2*sqrt(5));

m2 = 4/(1-sqrt(5));
c2 = (sqrt(5)-3)/(2*sqrt(5));

% E  = .25*2 + .125*( m1^3 *c1 * (()/((m1 -1 )*(1-m1))) + m2^3 *c2 * (1/((m2 -1 )*(1-m2)))^2);


E2 = .25*2 + .25*( m1^1 *c1 * (1/((m1 -1 )*(1-m1))) + m2^1 *c2 * (1/((m2 -1 )*(1-m2))));





