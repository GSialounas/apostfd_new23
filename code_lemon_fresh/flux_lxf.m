function f_h =flux_lxf(dist_x_pl,dt,uold)
F_pl = ( dist_x_pl./(2*dt).*(uold-circshift(uold,-1)) +.5*.5*(uold.^2 +circshift(uold,-1).^2  ));
F_min = ( circshift(dist_x_pl,1)./(2*dt).*(circshift(uold,1)-circshift(uold,0)) +.5*.5*(circshift(uold,1).^2 +circshift(uold,0).^2  ));
f_h = (1./dist_x_pl) .*(F_pl-F_min);
end