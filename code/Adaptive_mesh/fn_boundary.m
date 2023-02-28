function [h,v] = fn_boundary(h,v)

h(1) = h(2);
v(1) = v(2);
h(end) = h(end-1);
v(end) = v(end-1);