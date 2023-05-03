function f_out = fn_hat_periodic_elliptic(x)

f_out =zeros(size(x));
% f_out(x<.5) =  .25-abs(x(x<.5)-.25);
% f_out(x>=.5) = -.25+abs(x(x>=.5)-.75);

f_out(x<=1) = x-1;
x75 = x(x<.75);
f_out(x<.75) = .5 - x75;
x25 = x(x<.25);
f_out(x<.25) = x25;

end