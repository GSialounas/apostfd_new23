function u = fn_exact_step_burger(x,t,ul,ur)
u = zeros(size(x));

if ul>=ur
u(x<=t*(ul+ur)/2)=ul;
u(x>t*(ul+ur)/2) = ur;
else
u(x<=t*ul)=ul;
u(x>ul*t & x<=ur*t) = x(x>ul*t & x<=ur*t)/t;
u(x>t*ur) = ur;
end
end