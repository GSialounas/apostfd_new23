clc;
clear;
close all;



el1 = element(0,1,1,1e-6);
solution = [1,2];
bound = [2,3];
el1 = el1.assign_vals(solution,bound);

el1 = el1.refine();
el1.leftChild = el1.leftChild.refine();

% el1 = el1.assign_vals(solution,bound);

x= linspace(0,1,1001);
tol = 1e-6;
% mesh{1}.root.initTriangulation = 1;
tree_list = [];
root_ids = [];
for i =1:length(x)-1
    
    if i==1
     tree_list = [ Tree(element(x(i), x(i+1),i,tol),i,tol)];
     tree_list(i) = tree_list(i).assign_vals([sin(x(i)), sin(x(i+1))], exp(-100*(x(i)-.5).^2));
     root_ids(i) = i;
    else
        tree_list(end+1) = Tree(element(x(i), x(i+1),i,tol),i,tol);
        tree_list(i) = tree_list(i).assign_vals([sin(x(i)), sin(x(i+1))], exp(-100*(x(i)-.5).^2));
        root_ids(i) = i;
    end
   
end
tic;
% 
 M = mesh(tree_list,root_ids);
    t1 = toc;
 M.updateMesh();
 t2= toc;
M.refineGlobal();
t3=toc;

M.refineGlobal();
t4= toc;
M.refineGlobal();
t5= toc;
M.mark_ref_coar();
t6= toc;

[grid_new, bound_new] = M.getNewGridVals2();
t7 = toc;
t_arr = [t1, t2, t3, t4, t5, t6, t7];
t_arr = [t1, t_arr(2:end)-t_arr(1:end-1)];
M.ref_coarsen_local();
[grid_new, bound_new] = M.getNewGridVals2();
grid_spacing = grid_new(2:end)-grid_new(1:end-1);

figure;
plot(grid_new(1:end-1), grid_spacing,'b.');
M.mark_ref_coar();
M.ref_coarsen_local();
[grid_new, bound_new] = M.getNewGridVals2();
grid_spacing = grid_new(2:end)-grid_new(1:end-1);

figure;
plot(grid_new(1:end-1), grid_spacing,'b.');



% [grid_new, bound_new, ref_new] = M.getNewGridVals();

% bound_new = M.getBound();
% ref_new = M.getRefInd();
%  M.refineGlobal();
%  t1=toc;
% M.coarsenGlobal();
% M.coarsenGlobal();
% M.coarsenGlobal();
% M.coarsenGlobal();
% t2 = toc;

% t_arr = [t1, t2-t1];