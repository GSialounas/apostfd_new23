function [nonu_ders, nonu_inds] =  fn_nonu_non_periodic(x,dist_x_pl,y_i)


y_ext = [ y_i(2) y_i(1), y_i, y_i(end), y_i(end-1), y_i(end-2)];
y_im  = circshift(y_i,[0 1]);
y_ip  = circshift(y_i,[0 -1]);
y_ipp = circshift(y_i,[0,-2]);

h_im = x - circshift(x,1);
h_im(1) = h_im(2);
h_i  = circshift(x,-1) - x;
h_i(end) = h_i(end-1);
h_ip = circshift(x,-2) - circshift(x,-1);
h_ip(end-1) =h_ip(end);
H = circshift(x,-2)  - circshift(x,1);
H(1) = H(2);
H(end) = H(end-2);
H(end-1) = H(end-2);

% xi_m  =  x - dist_x_pl;
% xi    =  x;
% xi_p  =  x + dist_x_pl;
% xi_pp =  x + 2*dist_x_pl;

y_im_d = - ( ( (2*h_im + h_i).* H + h_im.* (h_im +h_i) )./( h_im.*(h_im + h_i).* H ) ).* y_im ...
       + (((h_im +h_i) .* H)./(h_im .* h_i.*(h_i +h_ip))) .* y_i...
       - ((h_im .* H)./((h_im + h_i).*h_i.*h_ip)) .* y_ip...
       + ((h_im.*(h_im +h_i))./((h_i +h_ip).*h_ip.*H)) .* y_ipp;

y_i_d = - ((h_i.*(h_i + h_ip))./(h_im .*(h_im + h_i).*H)).* y_im...
        + ((h_i.*(h_i + h_ip) - h_im.*(2*h_i +h_ip))./(h_im.*h_i.*(h_i + h_ip)) ).* y_i...
        + ((h_im.*(h_i +h_ip))./((h_im + h_i).*h_i.*h_ip)).* y_ip...
        - ((h_im.*h_i)./((h_i + h_ip).*h_ip.*H)).*y_ipp;
y_ip_d = ((h_i.*h_ip)./(h_im.*(h_im+h_i).*H)) .* y_im...
         - (h_ip.*(h_im+h_i)./(h_im.*h_i.*(h_i +h_ip))) .* y_i ...
         + (((h_im+ 2*h_i).*h_ip - (h_im + h_i).*h_i)./((h_im+h_i).*h_i.*h_ip)).*y_ip...
         + (((h_im + h_i).*h_i)./((h_i+h_ip).*h_ip.*H)).*y_ipp;
y_ipp_d = -(h_i+h_ip).*h_ip./(h_im.*(h_im+h_i).*H) .*y_im...
          + ((h_ip.*H)./(h_im.*h_i.*(h_i+h_ip))).*y_i...
          - ((h_i + h_ip).*H./((h_im +h_i).*h_i.*h_ip)).*y_ip...
          + (((2*h_ip +h_i).*H +h_ip.*(h_i +h_ip))./((h_i+h_ip).*h_ip.*H)).*y_ipp;
      
      
nonu_ders  = [y_im_d; y_i_d;y_ip_d;y_ipp_d];

B0 = (h_i + h_ip).^2 .* (abs(y_ip_d - y_i_d)./h_i - abs(y_i_d -y_im_d)./h_im).^2;
B1 = (h_im + h_i).^2 .*(abs(y_ipp_d- y_ip_d)./h_ip - abs(y_ip_d - y_ip)./h_i).^2;
end