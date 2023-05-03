function [s,ss] = fn_sum(n)


k = 3:n;

r1 =((1+sqrt(5))/4);
r2 = ((1-sqrt(5))/4);
c1 = (sqrt(5)+3)/(2*sqrt(5));
c2 = (sqrt(5)-3)/(2*sqrt(5));

sum1 = .125 .*(sqrt(5)+3)/(2*sqrt(5))*k .* ((1+sqrt(5))/4).^(k-3) ;
sum2 =.125 .*(sqrt(5)-3)/(2*sqrt(5))*k .* ((1-sqrt(5))/4).^(k-3);


s = .5+ sum(sum1) +sum(sum2);


sum1s = 1/(1-r1) * c1 *r1^(-3) * ((r1^3)/(1-r1)+2*r1^3);
sum2s = 1/(1-r2) * c2 *r2^(-3) * ((r2^3)/(1-r2)+2*r2^3);

ss = .5+1/8*(sum1s+sum2s);
end