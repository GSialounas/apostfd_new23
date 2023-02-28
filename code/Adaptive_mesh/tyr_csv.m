clc;
clear;
close all;
M = dlmread('../../paper/data/My_file.csv', ',', 2, 0);
M2= dlmread('../../paper/data/error_file_adapt.csv', ',', 2, 0);
figure
semilogy(M(:,1),M(:,2),'b',M2(:,1),M2(:,2),'r')
xlim([-.1, .5])

% data = load('../../paper/data/My_file.csv')
% x=  [1:1000]';
% y1=x.^2;
% y2 = x.^3;
% ss= ["time","error1","error2"];
% y= [x, y1,y2];
% 
% 
% 
% fid  = fopen('tablefile_gs.csv','w');
% fprintf(fid,'%s, %s, %s\n',ss(1), ss(2), ss(3))
% 
% for i =1:length(y)
% fprintf(fid, '%5.4f,%5.4f,%5.4f\n', y(i,:));
% end
% fclose(fid);