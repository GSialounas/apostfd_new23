function [t]  = assemble_and_solve_1d_residual(f_fn,S,u_fd)


global mesh uh fh wq nq xq

S(1,:)=0;
S(end,:)=0;
S(1,1) = 1;
S(end,end)=1;
n_vertices = mesh.n_vertices;
n_elem = mesh.n_elem;

fh = zeros(n_vertices,1);
tic

for el = 1:n_elem
    v_elem = mesh.elem_vertices(el,:);
    v1 = mesh.vertex_coordinates(v_elem(1));
    v2 = mesh.vertex_coordinates(v_elem(2));
    
    h = v2-v1;
    m = (v1+v2)/2;
    x_eval = m;
    
    integral_el = 0;
    for iq = 1:nq
      xiq = 0.5*(v2-v1)*xq(iq) + m;
      f_iq = f_fn(xiq);
      [IU_iq, IU_x_iq,IU_xx_iq] = fn_hermite_rec_ele(xiq,el,f_fn,u_fd);
      FIU_iq = 0.5*IU_iq.^2;
      
      
      % basis functions and their gradients
      phi_1_iq = (v2 - xiq)/(v2-v1);
      phi_1_x_iq = -1/(v2-v1);
      
      phi_2_iq = (xiq - v1)/(v2-v1);
      phi_2_x_iq = 1/(v2-v1);
      
      
      f_el = [wq(iq)*h*f_iq*phi_1_iq ; wq(iq)*h*f_iq*phi_2_iq];
      
      grad_el  = [wq(iq)*h*IU_x_iq *phi_1_x_iq ;  wq(iq)*h*IU_x_iq *phi_2_x_iq  ];
      
      
      % the lines that I need here are the following 
      f_el = [wq(iq)*h*IU_t_iq*phi_1_iq ; wq(iq)*h*IU_t_iq*phi_2_iq];
      
      grad_el  = [wq(iq)*h*FIU_iq *phi_1_x_iq ;  wq(iq)*h*FIU_iq *phi_2_x_iq  ];
      
      fh(v_elem) = fh(v_elem)+ f_el + grad_el;

    end
end
fh(1) =0;
fh(end) = 0;
uh=S\fh;

t=toc;
end