clc;
clear;
close all;


pdf_vec = zeros(1, 301);
p = 2/3;
q = 1/3;
for i = 0:300
    pdf_vec(i+1) = nchoosek(300,i)*p^i *q^(300-i);
end


x=  0:length(pdf_vec)-1;

figure
plot(x,pdf_vec,'bo')

cdf_vec = cumsum(pdf_vec);

figure
plot(x,cdf_vec,'ro')
% 
% t = 1:1:1000;
% figure
% for i =1:length(t)
%     plot(t(i),t(i),'b.')
%     hold on;
%     pause(0.0001)
% end
%     




x = linspace(0,32*pi,1001);

y= sin(1/16*pi*x);
ylog = exp(x);
y2 = cos(1/16*pi*x);
figure
yyaxis left
plot(x,y,'b');

yyaxis right
semilogy(x,ylog,'r')