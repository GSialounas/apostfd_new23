classdef mesh <handle
    properties 
        tree_list
        root_ids
        elem_list% = M.getElements();
        elem_list_flat %= M.getElements_flast();
    end
    
    methods
        function obj = mesh(tree_list,root_ids)
            obj.tree_list = tree_list;
            obj.root_ids = root_ids;
        end
        
        function obj= assign_rootID(obj)
            N_roots = length(mesh.tree_list);
            for i =1:N_roots
                
            end
        end
        

        
        
        function obj = ref_coarsen(obj)
            mesh.element_listref_coarsen()
        end
        
        function list  = getLeaves(obj)
            list = [];
            for i = 1:length(obj.tree_list)
                list= [list, obj.tree_list(i).getLeaves()];
            end
        end
        
        function list = get_element_vert(obj)
            list =[];
            for i = 1:length(obj.tree_list)
                x= obj.tree_list(i).root.getElement_vert([]);
                list = [list, x];
            end
        end
        
        function list = getElements(obj)
            list = cell(1, length(obj.tree_list));
%               list =getElements_arr(obj.tree_list);
            for i = 1:length(obj.tree_list)
%                 if  i == 1
                    list(i) = {obj.tree_list(i).getElements()};
%                 else
%                     list{end+1}=obj.tree_list(i).getElements();
%                 end
            end
        end
        function list = getElements_flat(obj)
            list = cell(1, length(obj.tree_list));
            for i = 1:length(obj.tree_list)
                if  i == 1
                    list{i} = obj.tree_list(i).getElements();
                else
                    list{i}= obj.tree_list(i).getElements();
                end
            end
            list = horzcat(list{:});
            
        end
        
        function obj = getElements_flat_time(obj)
%             tic;
            list = cell(1, length(obj.tree_list));
%             t1 = toc;

            for i = 1:length(obj.tree_list)
                list{i} = getElements(obj.tree_list(i));
            end
%             t2 = toc;
            list = horzcat(list{:});
%             t3 = toc;

%             t_arr =[t1, t2-t1, t3-t2];
            
        end

        function t_arr = getElements_flat_time_arr(obj)
            tic;
%             list = cell(1, length(obj.tree_list));
            t1 = toc;
            list = getElements_arr(obj.tree_list);
%             for i = 1:length(obj.tree_list)
%                 list{i} = getElements(obj.tree_list(i));
%             end
            t2 = toc;
            list = horzcat(list{:});
            t3 = toc;

            t_arr =[t1, t2-t1, t3-t2];
            
        end

%         function C = getElements_flat_d(tree_list)
%             C.elem_list_flat{i}
%         end
%         function list = getElements_flat(obj)
%             list = {};
%             for i = 1:length(obj.tree_list)
%                 if  i == 1
%                     list = obj.tree_list(i).getElements();
%                 else
%                     list=[list, obj.tree_list(i).getElements()];
%                 end
%             end
%             
%         end
        function obj = updateMesh(obj)
%             tic;
            obj.elem_list = getElements_arr(obj.tree_list);%obj.getElements();

%             t1 = toc;
            obj.elem_list_flat = horzcat(obj.elem_list{:});

%             t2 = toc;
%             t_arr = [t1, t2-t1];
        end
        
        function obj = asgnVals(obj,new_vals_soln, new_vals_bound)
%             list_ele = obj.getElements_flat();
            % this is just for periodic boundary conds
%             if length(new_soln)~= (length(obj.elem_list_flat)+1)
%                 disp("do you want periodic bcs or not.  Check because length ofnew vals doesn't match")
%             else
                

            for i = 1:length(obj.elem_list_flat)
                if i < length(obj.elem_list_flat)
                    obj.elem_list_flat(i).solution =[new_vals_soln(i),new_vals_soln(i+1)];
                    obj.elem_list_flat(i).bound = new_vals_bound(i);%[new_vals_bound(i),new_vals_bound(i+1)];
                else
                    obj.elem_list_flat(i).solution =[new_vals_soln(i),new_vals_soln(1)];
                    obj.elem_list_flat(i).bound = new_vals_bound(i);%[new_vals_bound(i),new_vals_bound(1)];
                end

            end
        end
         
        
        function grid = getGrid(obj)
            % creates grids from leaf elements
            grid = zeros(1,length(obj.elem_list_flat));
            for i = 1:length(grid)
                grid(i)= obj.elem_list_flat(i).leftEnd;
            end
        end
        
        function bound_arr = getBound(obj)
            % creates grids from leaf elements
            bound_arr = zeros(1,length(obj.elem_list_flat));
            for i = 1:length(bound_arr)
                bound_arr(i)= obj.elem_list_flat(i).bound;
            end
        end
        
        function ref_ind = getRefInd(obj)
            % creates grids from leaf elements
            ref_ind = zeros(1,length(obj.elem_list_flat));
            for i = 1:length(ref_ind)
                ref_ind(i)= obj.elem_list_flat(i).ref_indicator;
            end
        end
        
        function [grid, bound_arr, ref_ind] = getNewGridVals(obj)
            grid = zeros(1,length(obj.elem_list_flat));
            bound_arr = zeros(1,length(obj.elem_list_flat));
            ref_ind = zeros(1,length(obj.elem_list_flat));
            for i = 1:length(grid)
                grid(i)= obj.elem_list_flat(i).leftEnd;
                bound_arr(i)= obj.elem_list_flat(i).bound;
                ref_ind(i)= obj.elem_list_flat(i).ref_indicator;

            end

        end
        
        function [grid, bound_arr] = getNewGridVals2(obj)
            grid = zeros(1,length(obj.elem_list_flat));
            bound_arr = zeros(1,length(obj.elem_list_flat));
            for i = 1:length(grid)
                grid(i)= obj.elem_list_flat(i).leftEnd;
                bound_arr(i)= obj.elem_list_flat(i).bound;

            end

        end
        
%         function bound_arr = getBound(obj)
%             % creates grids from leaf elements
%             bound_arr = zeros
%             for i = 1:length(obj.tree_list)
%                 grid= [grid, obj.tree_list(i).getLeaves()];
%             end
%         end
        
        
        function obj= refineGlobal(obj)
%             tic;
%             obj.updateMesh();
%             t1 = toc;
            ele_flat_length = length(obj.elem_list_flat);
%             t2=toc;
            
            for i = 1:ele_flat_length
%                 obj.elem_list_flat(i).refine();
                obj.elem_list_flat(i).refine();

            end
%             t3= toc;
            % update properties of the mesh
            obj.updateMesh();
%             t4 = toc;
            
%             t_arr = [t1, t2-t1, t3-t2, t4-t3];
            
        end
        
        function obj = coarsenGlobal(obj)
%             obj.updateMesh();
            ele_flat_length = length(obj.elem_list_flat);
            
            for i = 1:ele_flat_length
                
                obj.elem_list_flat(i).coarsen();
            end
            
            % update properties of the mesh
            obj.updateMesh();
        end
        
        function obj = mark_ref_coar(obj)
%             obj.updateMesh();
            ele_flat_length = length(obj.elem_list_flat);
            max_bound = 0;
            for i = 1:ele_flat_length
                max_bound = max(max_bound, obj.elem_list_flat(i).bound);
            end
            
            for i = 1:ele_flat_length
                if obj.elem_list_flat(i).bound > 0.9 * max_bound
                    obj.elem_list_flat(i).ref_indicator = 1;
                elseif obj.elem_list_flat(i).bound < 0.1 * max_bound
                    obj.elem_list_flat(i).ref_indicator = -1;
                else
                     obj.elem_list_flat(i).ref_indicator=0;
                end
            end
            
            obj.updateMesh();
            
            
        end
        
        
        
        function obj = ref_coarsen_local(obj)
            
%             obj.updateMesh();
            ele_flat_length = length(obj.elem_list_flat);
            
            for i = 1:ele_flat_length
                if (obj.elem_list_flat(i).ref_indicator ==1)
                    obj.elem_list_flat(i).refine();
                elseif obj.elem_list_flat(i).ref_indicator == -1
                    obj.elem_list_flat(i).coarsen();
                else
                end
            end
            
            % update properties of the mesh
            obj.updateMesh();
        end
%         
    end
    
end