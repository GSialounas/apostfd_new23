function F_out = flux_ww_dambreak(dist_x_pl,dt,u,F)
% This is not lax-wendroff but rather white wendroff.  It's a modified
% version of lax-wendforff which maintains second order convergence for
% irregular grids

u = [u(:,1), u ,u(:,end)]; %[ h; hv]
F = [F(:,1), F, F(:,end)]; % [F_h; F_hv]
dist_x_pl = [dist_x_pl(1), dist_x_pl, dist_x_pl(end)];

u_min = zeros(size(u));
u_pl = zeros(size(u));

h= u(1,:);
v= u(2,:)./u(1,:);
g = 9.81;

dx = dist_x_pl(1);
for j  = 2:length(u)-1
    u_pl(:,j-1)  = .5 * (u(:,j+1) + u(:,j)) - .5 * (dt/dx) * (F(:, j+1) - F(:,j));
    u_min(:,j-1) = .5 * (u(:,j) + u(:,j-1)) - .5 * (dt/dx) * (F(:,j) - F(:,j-1));
    
end

h_pl = u_pl(1,:);
v_pl = u_pl(2,:)./u_pl(1,:);

h_min = u_min(1,:);
v_min = u_min(2,:)./u_min(1,:);

[h_min, hv_min, F_h_min, F_hv_min] = dependent(h_min,v_min,g);
[h_pl, hv_pl, F_h_pl, F_hv_pl] = dependent(h_pl,v_pl,g);


% 
% for j = 2:length(u)-1
%     u_min(:,j-1) = (1/(dist_x_pl(j-1)+dist_x_pl(j))) * (dist_x_pl(j-1)*u(:,j) +dist_x_pl(j)*u(:,j-1))+...
%              dt/(dist_x_pl(j)+dist_x_pl(j-1)) * (F(:,j)-F(:,j-1));
%     u_pl(:,j-1)  = (1/(dist_x_pl(j)+dist_x_pl(j+1))) * (dist_x_pl(j)*u(:,j+1) + dist_x_pl(j+1)* u(:,j))+...
%              (dt/(dist_x_pl(j+1)+dist_x_pl(j))) * (F(:,j+1)-F(:,j));
%          
% end
% 

F_min = [F_h_min;F_hv_min];
F_pl = [F_h_pl; F_hv_pl];
F_out = (1/dx)*(F_pl(:,2:end-1)-F_min(:,2:end-1));


end