function u  = fn_ex_transport_stp(x,t)

u=zeros(size(x));
u(x<t)=1;
end
