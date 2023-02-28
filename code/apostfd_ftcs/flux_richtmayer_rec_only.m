function f_h = flux_richtmayer_rec_only(dist_x_pl,dt,u)
    u_pl  = .5 * (circshift(u,-1) + u) - (dt./(2*dist_x_pl)) .* ((circshift(u,-1)) - u);
    u_min = .5 * (circshift(u,1) + u) - (dt./(2*dist_x_pl)) .* (u - circshift(u,1));
    f_h = (1./dist_x_pl) .*  (u_pl  - u_min);
end