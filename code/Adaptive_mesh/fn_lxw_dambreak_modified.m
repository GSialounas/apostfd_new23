function [f_h_out, h_new] = fn_lxw_dambreak_modified(h,v,dist_x_pl, dt,g)

f_h_out = zeros(2, length(h));
h_new = zeros(2, length(h));
h = [h(1) h h(end)];
v = [v(1) v v(end)];
h_min = zeros(size(h));
v_min = zeros(size(v));
h_pl = zeros(size(h));
v_pl =  zeros(size(v));

[un1, un2, F1, F2] = dependent(h,v,g);
un1_new = zeros(size(un1));
un2_new  =zeros(size(un2));
un1_pl = zeros(size(un1));
un1_min = zeros(size(un1));
un2_pl = zeros(size(un1));
un2_min = zeros(size(un1));
dist_x_pl = [dist_x_pl(1) dist_x_pl dist_x_pl(end)];

for j = 2:length(h)-1
    un1_pl(j) = .5*(un1(j+1) + un1(j)) - .5 * dt/(dist_x_pl(j)) *(F1(j+1)- F1(j));
    un2_pl(j) = .5*(un2(j+1) + un2(j)) - .5 * dt/(dist_x_pl(j)) * (F2(j+1) - F2(j));
    
    un1_min(j) = .5 * (un1(j) + un1(j-1)) - .5 *dt/(dist_x_pl(j-1)) * (F1(j)- F1(j-1));
    un2_min(j) = .5 * (un2(j) + un2(j-1)) - .5 *dt/(dist_x_pl(j-1)) * (F2(j) - F2(j-1));
    
    
    
%     %%%%%%%%
%     un1_pl(j) =  (dist_x_pl(j-1)*un1(j+1)+ dist_x_pl(j)*un1(j))/(dist_x_pl(j)+dist_x_pl(j-1))...
%         - .5 * (dt/dist_x_pl(j)) * (F1(j+1) - F1(j));
%     un2_pl(j) = (dist_x_pl(j-1)*un2(j+1)+ dist_x_pl(j)*un2(j))/(dist_x_pl(j)+dist_x_pl(j-1))...
%         - .5 * (dt/dist_x_pl(j)) * (F2(j+1) - F2(j));
%     
%     if j==2
%         un1_min(j) = .5 * (un1(j)+un1(j-1)) - .5 * (dt/dist_x_pl(j-1)) * (F1(j) - F1(j-1));
%         un2_min(j) = .5 * (un2(j)+un2(j-1)) - .5 * (dt/dist_x_pl(j-1)) * (F2(j) - F2(j-1));
%         
%       
%     else
%         un1_min(j) = (dist_x_pl(j-2)*un1(j)+dist_x_pl(j-1)*un1(j-1))/(dist_x_pl(j-2)+dist_x_pl(j-1)) ...
%             - .5 * (dt/dist_x_pl(j-1)) * (F1(j) - F1(j-1));
%         un2_min(j) = (dist_x_pl(j-2)*un2(j)+dist_x_pl(j-1)*un2(j-1))/(dist_x_pl(j-2)+dist_x_pl(j-1))...
%             - .5 * (dt/dist_x_pl(j-1)) * (F2(j) - F2(j-1));
%     end
%     
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
    un1_new(j) = un1(j) - 2*(dt/(dist_x_pl(j)+dist_x_pl(j-1))) * (F1_pl(j) - F1_min(j));
    un2_new(j) = un2(j) - 2*(dt/(dist_x_pl(j)+dist_x_pl(j-1))) * (F2_pl(j) - F2_min(j));
    f_h_out(:, j-1) = [2/(dist_x_pl(j)+dist_x_pl(j-1)) * (F1_pl(j) - F1_min(j));2/(dist_x_pl(j)+dist_x_pl(j-1))  * (F2_pl(j) - F2_min(j))];
    
    %%%%%%
%     un1_new(j) = un1(j) - (dt/dist_x_pl(j-1)) * (F1_pl(j) - F1_min(j));
%     un2_new(j) = un2(j) - (dt/dist_x_pl(j-1)) * (F2_pl(j) - F2_min(j));
%     f_h_out(:, j-1) = [1/dist_x_pl(j-1) * (F1_pl(j) - F1_min(j));1/dist_x_pl(j-1) * (F2_pl(j) - F2_min(j))];
end

h = un1_new;
v= un2_new./un1_new;
[h,v] = fn_boundary(h,v);
[un1, un2, F1, F2] = dependent(h,v,g);

 f_h_out = [1./dist_x_pl; 1./dist_x_pl] .* [F1(2:end)- F1(1:end-1); F2(2:end) - F2(1:end-1)];

h_new = [un1(2:end-1); un2(2:end-1)];



end