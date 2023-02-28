function F_out =flux_lxf_dambreak(dist_x_pl,dt,uold,F)
% f_pl = ( dist_x_pl./(2*dt)).*(uold-circshift(uold,-1)) +.5*(F +circshift(F,-1));
% f_min = ( circshift(dist_x_pl,1)./(2*dt)).*(circshift(uold,1)-uold) +.5*(circshift(F,1) + F);
uold = [uold(1), uold ,uold(end)];
F = [F(1), F, F(end)];
dist_x_pl = [dist_x_pl(1), dist_x_pl, dist_x_pl(end)];
f_pl = ( dist_x_pl(2:end)./(2*dt)).*(uold(2:end-1)-uold(3:end)) +.5*(F(2:end-1) +F(3:end));
f_min = (dist_x_pl(1:end-1)./(2*dt)).*(uold(1:end-2)-uold(2:end-1)) +.5*(F(1:end-2) + F(2:end-1));

F_out = (1./dist_x_pl(1:end-1)) .*(f_pl-f_min);
end
