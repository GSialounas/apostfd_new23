function [nonu_ders, nonu_inds] =  fn_nonu_ders_nonu(x_v,y_i,L_domain)

y_im  = circshift(y_i,[0 1]);
y_ip  = circshift(y_i,[0 -1]);
y_ipp = circshift(y_i,[0,-2]);

% y_im = y_i(:,1:end-3);
% y_ipp = y_i(:,4:end);
% y_ip = y_i(:,3:end-1);
% y_i = y_i(:,2:end-2);


% y_im  = circshift(y_i,[0 1]);
% y_i  = circshift(y_i,[0 0]);
% y_ip = circshift(y_i,[0,-1]);
% y_ipp = circshift(y_i,[0 -2]);
% 
% xi_m = x_v(1:end-3);
% xi_p = x_v(3:end-1);
% xi = x_v(1:end-2);
% xi_pp = x_v(4:end);



% h_i = xi_p -xi;
% h_im = xi - xi_m;
% h_imm = xi_m -xi_mm;



h_i = circshift(x_v,-1)- x_v;
h_i(end)=L_domain-x_v(end);

h_im = circshift(h_i,1);
h_ip = circshift(h_i,-1);
% h_im = circshift(x,1) - circshift(x,2);
% h_im(2)=L_domain-x(end);

% h_ip  = circshift(x,-1) - x;
% h_ip(end) = h_ip(end-1);


% 
% h_ipp = circshift(x,-2) - circshift(x,-1);
% h_ipp(end-1) =h_ipp(end);

% h_ip = circshift(x,-1) - circshift(x,0);
% h_ip(end) =L_domain-x(end);%h_ip(1);
% 


H = h_im+h_i+h_ip;% circshift(x,-1)  - circshift(x,2);

% 
% H(1) = H(3);
% H(2) = H(3);
% H(end) =H(end-1);

% xi_m  =  x - dist_x_pl;
% xi    =  x;
% xi_p  =  x + dist_x_pl;
% xi_pp =  x + 2*dist_x_pl;

y_im_d = - ( ( (2*h_im + h_i).* H + h_im.* (h_im +h_i) )./( h_im.*(h_im + h_i).* H ) ).* y_im ...
       + (((h_im +h_i) .* H)./(h_im .* h_i.*(h_i +h_ip))) .* y_i...
       - ((h_im .* H)./((h_im + h_i).*h_i.*h_ip)) .* y_ip...
       + ((h_im.*(h_im +h_i))./((h_i +h_ip).*h_ip.*H)) .* y_ipp;

y_i_d = - ((h_i.*(h_i + h_ip))./(h_im .*(h_im + h_i).*H)).* y_im...
        + ((h_i.*(h_i + h_ip) - h_ip.*(2*h_i +h_ip))./(h_im.*h_i.*(h_i + h_ip)) ).* y_i...
        + ((h_im.*(h_i +h_ip))./((h_im + h_i).*h_i.*h_ip)).* y_ip...
        - ((h_im.*h_i)./((h_i + h_ip).*h_ip.*H)).*y_ipp;
    
y_ip_d = ((h_i.*h_ip)./(h_im.*(h_ip+h_i).*H)) .* y_im...
         - (h_ip.*(h_im+h_i)./(h_im.*h_i.*(h_i +h_ip))) .* y_i ...
         + (((h_im+ 2*h_i).*h_ip - (h_im + h_i).*h_i)./((h_im+h_i).*h_i.*h_ip)).*y_ip...
         + (((h_im + h_i).*h_i)./((h_i+h_ip).*h_ip.*H)).*y_ipp;
     
y_ipp_d = -(h_i+h_ip).*h_ip./(h_im.*(h_im+h_i).*H) .*y_im...
          + ((h_ip.*H)./(h_im.*h_i.*(h_i+h_ip))).*y_i...
          - ((h_i + h_ip).*H./((h_im +h_i).*h_i.*h_ip)).*y_ip...
          + (((2*h_ip +h_i).*H +h_ip.*(h_i +h_ip))./((h_i+h_ip).*h_ip.*H)).*y_ipp;
      
      
nonu_ders  = [y_im_d; y_i_d;y_ip_d;y_ipp_d];

B0 = (h_i + h_ip).^2 .* (abs(y_ip_d - y_i_d)./h_i - abs(y_i_d -y_im_d)./h_im).^2;
B1 = (h_im + h_i).^2 .*(abs(y_ipp_d- y_ip_d)./h_ip - abs(y_ip_d - y_i_d)./h_i).^2;
nonu_inds = [B0;B1];
end