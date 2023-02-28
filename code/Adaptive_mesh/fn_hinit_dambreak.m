function h_initial = fn_hinit_dambreak(x)

h_initial= zeros(1,length(x));
x_dam = 30;
h_initial(x>x_dam)= .1;
h_initial(x<=x_dam) = .2;
end