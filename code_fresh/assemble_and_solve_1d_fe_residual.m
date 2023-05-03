function [R_i,L1Rt,L1Rt_arr]  = assemble_and_solve_1d_fe_residual(S,x,uold,u,evalt,tj,dt,ul,ur)


global mesh uh fh wq nq xq

L1Rt = 0;
L1Rt_arr = zeros(size(x));

h = x(2)-x(1);

S(1,:)=0;
S(end,:)=0;
S(1,1) = 1;
S(end,end)=1;
n_vertices = mesh.n_vertices;
n_elem = mesh.n_elem;

fh = zeros(n_vertices,1); % n_vertices -1 because of periodic domain
tic

% temporal reconstruction
c0 = uold ;
c1 = (1/dt) *(u-uold);

diff_t =evalt-tj;

Ut = c0 + c1*diff_t;
Ut_t = c1;

% Spatial reconstruction
c_0_ts = Ut; 
c_1_ts = 1./h.*(circshift(Ut,-1)-circshift(Ut,0));


c_0_ts_t = Ut_t; 
c_1_ts_t = 1./h.*(circshift(Ut_t,-1)-circshift(Ut_t,0));

for el = 1:n_elem
    v_elem = mesh.elem_vertices(el,:);
    v1 = mesh.vertex_coordinates(v_elem(1));
    v2 = mesh.vertex_coordinates(v_elem(2));
    
    h = v2-v1;
    m = (v1+v2)/2;
    x_eval = m;
    
    integral_el = 0;
    
    % temporal reconstruction 
    c_0_ts_elem  = c_0_ts(v_elem(1));
    c_1_ts_elem = c_1_ts(v_elem(1));
    c_0_ts_t_elem = c_0_ts_t(v_elem(1));
    c_1_ts_t_elem = c_1_ts_t(v_elem(1));

    
    for iq = 1:nq
      xiq = 0.5*(v2-v1)*xq(iq) + m;
      diff_x = xiq-x(v_elem(1));
      
      
%       f_iq = f_fn(xiq);
%      
%       [IU_iq, IU_x_iq,IU_xx_iq] = fn_hermite_rec_ele(xiq,el,f_fn,u_fd);
%       FIU_iq = 0.5*IU_iq.^2;
      
      IU = c_0_ts_elem + c_1_ts_elem .* diff_x ;
%       FIU = 0.5 *IU.^2;
      FIU = IU;

      FIU_x = 2*(c_0_ts_elem+c_1_ts_elem.*diff_x).*(c_1_ts_elem);
      IU_t = c_0_ts_t_elem + c_1_ts_t_elem.*diff_x;
      
      
      % basis functions and their gradients
      phi_1_iq = (v2 - xiq)/(v2-v1);
      phi_1_x_iq = -1/(v2-v1);
      
      phi_2_iq = (xiq - v1)/(v2-v1);
      phi_2_x_iq = 1/(v2-v1);

      
      % the lines that I need here are the following 
      f_el = [wq(iq)*h*IU_t*phi_1_iq ; wq(iq)*h*IU_t*phi_2_iq];
      
      grad_el  = [wq(iq)*h*FIU *phi_1_x_iq ;  wq(iq)*h*FIU *phi_2_x_iq  ];
      
      fh(v_elem) = fh(v_elem)+ f_el - grad_el;

    end
end
% fh(1) =ul;
% fh(end) = ur;
uh=S\fh(1:end-1);

% Now we have R_i =R(x_i) from the uh.  Let us now obtain ||R||_L1(W)
R_i = uh;
c_0_R = R_i;
c_1_R = (1/h)*(circshift(R_i,-1)-R_i);

mask = ones(size(x));
mask(abs(x)>1)=0;
for iq = 1:nq
    xiq = .5*h*xq(iq) + x +h/2;
    diff_x = xiq-x;
    
    IR = c_0_R' + c_1_R' .* diff_x ;
    
    
    L1Rt_arr = L1Rt_arr + (wq(iq)*h*mask.*abs(IR));
    L1Rt = L1Rt + sum(wq(iq)*h*mask.*abs(IR));

end

t=toc;
end