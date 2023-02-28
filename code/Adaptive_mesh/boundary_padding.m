function [h,v] = boundary_padding(h,v)
% we pad the functions to create the fictitious node necessary for neumann
% bcs with zero gradient
h = [h(1) h h(end)];
v = [v(1) v v(end)];
end