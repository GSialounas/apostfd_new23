function [nonu_ders, nonu_inds] =  fn_nonu_ders_nonu_non_periodic(x_v,y_i,L_domain)


%This are the appropriate non-periodic smoothness indicators
% y_im  = circshift(y_i,[0 1]);
% y_ip  = circshift(y_i,[0 -1]);
% y_ipp = circshift(y_i,[0,-2]);
% 
% y_imm  = circshift(y_i,[0 2]);
% y_i  = circshift(y_i,[0 0]);
% y_ip = circshift(y_i,[0,-1]);
% y_im = circshift(y_i,[0 1]);
% 
% 
% 
% 
% % h_imm = circshift(x,1) - circshift(x,2);
% % h_imm(2)=L_domain-x(end);
% % 
% % 
% % 
% % h_im  = circshift(x,0) - circshift(x,1);
% % h_im(1) = L_domain-x(end);
% % % 
% 
% 
% h_i = circshift(x,-1) - circshift(x,0);
% h_i(end) =L_domain-x(end);%h_i(1);
% 
% h_im = circshift(h_i,1);
% h_imm = circshift(h_i,2);
y_imm = y_i(:,1:end-3);
y_im = y_i(:,2:end-2);
y_ip = y_i(:,4:end);
y_i = y_i(:,3:end-1);

xi_mm = x_v(1:end-3);
xi_m = x_v(2:end-2);
xi_p = x_v(4:end);
xi = x_v(3:end-1);



h_i = xi_p -xi;
h_im = xi - xi_m;
h_imm = xi_m -xi_mm;


H = h_imm+h_im+h_i;% circshift(x,-1)  - circshift(x,2);

% 
% H(1) = H(3);
% H(2) = H(3);
% H(end) =H(end-1);

% xi_m  =  x - dist_x_pl;
% xi    =  x;
% xi_p  =  x + dist_x_pl;
% xi_pp =  x + 2*dist_x_pl;

y_imm_d = - ( ( (2*h_imm + h_im).* H + h_imm.* (h_imm +h_im) )./( h_imm.*(h_imm + h_im).* H ) ).* y_imm ...
       + (((h_imm +h_im) .* H)./(h_imm .* h_im.*(h_im +h_i))) .* y_im...
       - ((h_imm .* H)./((h_imm + h_im).*h_im.*h_i)) .* y_i...
       + ((h_imm.*(h_imm +h_im))./((h_im +h_i).*h_i.*H)) .* y_ip;

y_im_d = - ((h_im.*(h_im + h_i))./(h_imm .*(h_imm + h_im).*H)).* y_imm...
        + ((h_im.*(h_im + h_i) - h_i.*(2*h_im +h_i))./(h_imm.*h_im.*(h_im + h_i)) ).* y_im...
        + ((h_imm.*(h_im +h_i))./((h_imm + h_im).*h_im.*h_i)).* y_i...
        - ((h_imm.*h_im)./((h_im + h_i).*h_i.*H)).*y_ip;
    
y_i_d = ((h_im.*h_i)./(h_imm.*(h_i+h_im).*H)) .* y_imm...
         - (h_i.*(h_imm+h_im)./(h_imm.*h_im.*(h_im +h_i))) .* y_im ...
         + (((h_imm+ 2*h_im).*h_i - (h_imm + h_im).*h_im)./((h_imm+h_im).*h_im.*h_i)).*y_i...
         + (((h_imm + h_im).*h_im)./((h_im+h_i).*h_i.*H)).*y_ip;
     
y_ip_d = -(h_im+h_i).*h_i./(h_imm.*(h_imm+h_im).*H) .*y_imm...
          + ((h_i.*H)./(h_imm.*h_im.*(h_im+h_i))).*y_im...
          - ((h_im + h_i).*H./((h_imm +h_im).*h_im.*h_i)).*y_i...
          + (((2*h_i +h_im).*H +h_i.*(h_im +h_i))./((h_im+h_i).*h_i.*H)).*y_ip;
      
      
nonu_ders  = [y_imm_d; y_im_d;y_i_d;y_ip_d];

B0 = (h_im + h_i).^2 .* (abs(y_i_d - y_im_d)./h_im - abs(y_im_d -y_imm_d)./h_imm).^2;
B1 = (h_imm + h_im).^2 .*(abs(y_ip_d- y_i_d)./h_i - abs(y_i_d - y_im_d)./h_im).^2;
nonu_inds = [B0;B1];
end