function u = fn_exact_step_combo_burger(x,t,ul,ur, xl, xr)
u = zeros(size(x));


u(x<=xl + t*(ul)/2)=ul;
u(x>(ul*t/2+xl) & x<=(ul*t+xr)) = 0;
u(x>(xr) & x<=(xr+ur*t)) = (x(x>(xr) & x<=(xr+ur*t))-xr)/t;
u(x>(xr+t*ur)) = ur;

end