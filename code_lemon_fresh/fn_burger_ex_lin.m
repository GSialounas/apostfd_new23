function [u,u_x,u_t] = fn_burger_ex_lin(x,t)


u_x=  zeros(size(x));
u_t = zeros(size(x));
if t<=1
    u = zeros(size(x));
    
    mask_1 = x>-1 & x<t;%(-1<(x-t)/(1+t)) & ((x-t)/(1+t)<0);
    %
    u(mask_1) = (1+x(mask_1))./(1+t);
    u_x(mask_1) = (1)./(1+t);
    u_t(mask_1) = -(1+x(mask_1))./((1+t).^2);
    
    mask_2 = x>=t;%(((x-t)/(1-t))>0) & (((x-t)/(1-t))<1);
    
    u(mask_2) = (1-x(mask_2))./(1-t);
    u_x(mask_2) = (-1)./(1-t);
    u_t(mask_2) = (1-x(mask_2))./((1-t).^2);
    
    u(x>1) = 0;
    u_x(x>1) = 0;
    u_t(x>1) = 0;

else
    u = zeros(size(x));
    
    mask_1 = (x>-1) & (x<=(sqrt(2)*sqrt(1+t)-1));%(0.5*t+0.5);%(-1<(x-t)/(1+t)) & ((x-t)/(1+t)<0);
    %
    u(mask_1) = (1+x(mask_1))./(1+t);
    u_x(mask_1) = (1)./(1+t);
    u_t(mask_1) = -(1+x(mask_1))./((1+t).^2);

    mask_2 = (x>(sqrt(2)*sqrt(1+t)-1));%0.5*t+0.5;%(((x-t)/(1-t))>0) & (((x-t)/(1-t))<1);
    
    u(mask_2) = 0;
    u_x(mask_2) = 0;
    u_t(mask_2) = 0;
end


end