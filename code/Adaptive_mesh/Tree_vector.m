classdef Tree_vector <handle
      
    
    properties
        root
        leaves
        root_id
        tol
    end
    
    methods
        function obj = Tree_vector(root,root_id,tol)
            obj.root = root;
            obj.root_id = root_id;
            obj.tol = tol;
        end
            
        function obj = refine(obj)
            obj.getLeaves.refine()            
        end
        
        function obj = assign_vals(obj, solution, bound)
            obj.root = obj.root.assign_vals(solution, bound);
        end
        
        function obj = insert_indicator(obj, data_error, data_indicator)
            % assign value of indicator for leaf elements
            % also assigns value of error
        end
        
        function obj = coarsen(obj)
            obj.getLeaves.coarsen()
        end
        
        function list_of_leaves = getLeaves(obj)
            list_of_leaves = unique(obj.root.getLeaves([]));
            
        end
        
        function list_of_elements = getElements(obj)
            list_of_elements  = obj.root.getElements([]);

        end
        function list_of_elements_arr = getElements_arr(obj)
            list_of_elements_arr = cell(1, numel(obj));
            for i = 1:numel(obj)
                list_of_elements_arr{i}  = obj(i).root.getElements([]);
            end

        end
        
        function list_of_elements_vert = getElement_vert(obj)
            list_of_elements_vert = obj.root.getElement_vert({});
        end
        
        function list_of_element_soln = getElement_soln(obj)
            list_of_element_soln = obj.root.getElement_soln({});
        end
        
        function list_of_element_bound = getElement_bound(obj)
            list_of_element_bound = obj.root.getElement_bound({});
        end
              
    end
    
end