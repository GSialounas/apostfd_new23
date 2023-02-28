classdef element < handle
    properties
        ele_ID
        
        bound
        solution
        
        tol
        macro_element

        leftChild
        rightChild
        
        leftEnd
        rightEnd
        leftRight_ind
        
        % ref_indicator:
        %  -1 : coarsen
        %  0  : do nothing
        %  +1 : refine
        ref_indicator
        
        
        parent
        initTriangulation=1
        refLevel = 0
        leaf
        
    end
    
    methods
        function obj = element(leftEnd, rightEnd,macro_element,tol)           
            obj.leftEnd = leftEnd;
            obj.rightEnd= rightEnd;
            obj.macro_element = macro_element;
            obj.tol= tol;
        end
        
        function obj = refine(obj)
            % grid data
            
            x_l = obj.leftEnd;
            x_m = .5*(obj.leftEnd + obj.rightEnd);
            x_r  = obj.rightEnd;
            
            
            tol = obj.tol;
            
            obj1 = element(x_l, x_m, obj.macro_element,tol);
            obj1.parent = obj;
            obj1.macro_element = obj1.parent.macro_element;
            obj1.initTriangulation = 0;
            obj1.refLevel = obj1.parent.refLevel +1;
            
            obj1.solution = [obj1.parent.solution(1), ...
                             .5*(obj1.parent.solution(1)+obj1.parent.solution(2))];
            %obj1.bound = [obj1.parent.bound(1), .5*(obj1.parent.bound(1)+obj1.parent.bound(2))];
            obj1.bound = [];% exp(-100*(x_l-.5).^2);
            obj1.leftRight_ind = "L";
            
            obj2 = element(x_m, x_r,obj.macro_element,tol);          
            obj2.parent = obj;
            obj2.macro_element = obj2.parent.macro_element;
            obj2.initTriangulation = 0;
            obj2.refLevel = obj2.parent.refLevel +1;
            obj2.solution =   [.5*(obj2.parent.solution(1)+obj2.parent.solution(2)),...
                                   obj2.parent.solution(2)];
            obj2.bound = [];%exp(-100*(x_m-.5).^2);
            obj2.leftRight_ind= "R";
            %obj2.bound =   [.5*(obj2.parent.bound(1)+obj2.parent.bound(2)), obj2.parent.bound(2)];

            
            obj.leftChild = obj1;
            obj.rightChild = obj2;
            
        end
         
        function  obj = coarsen(obj)
%             if isempty(obj.parent.leftChild)
%                 disp('leftChild is empty for element with macro')
%                 obj.parent.macro_element
%             elseif isempty(obj.parent.rightChild)
%                 
%                 disp('right is empty for element with macro')
%                 obj.parent.macro_element
%             else
%                 disp('nothing of the two for element with macro')
%                 obj.parent.macro_element
% %             end
%             if obj.leftRight_ind == 'L'
%                 % the iterator will reach this cell before its sibling
%                 obj.parent.solution = [obj.solution(1), obj. parent.rightChild.solution(2)];
% %                 disp('on left at some point')
%                 obj.parent.leftChild = [];
%                 obj.parent.rightChild = [];
%             else
%                 % in this case 
%                 obj.parent.leftChild.solution;
%                 obj.parent.rightChild.solution;
%                 obj.parent.solution = [obj.parent.leftChild.solution(1), obj.parent.rightChild.solution(2)];
                %                 disp('on left at some point')
                obj.parent.leftChild = [];
                obj.parent.rightChild = [];            %    disp('on right at some point')
%             end
%             obj.parent.solution = [obj.parent.leftChild.solution(1), obj.parent.rightChild.solution(2)];
            
            
        end
        
        function list = getLeaves(obj,list)
            if (isempty(obj.leftChild) && isempty(obj.rightChild))
%                 list= [list, obj.leftEnd, obj.rightEnd];
                list= [list, obj.leftEnd];

            else
                list = getLeaves(obj.leftChild,list);
                list = getLeaves(obj.rightChild,list);
            end
        end
        
        function list = getElement_vert(obj,list)
            if (isempty(obj.leftChild) && isempty(obj.rightChild))
                 %list{end+1}= [obj.leftEnd, obj.rightEnd];
                 list = [list, obj.leftEnd];

                
            else
                list = getElement_vert(obj.leftChild,list);
                list = getElement_vert(obj.rightChild,list);
            end
        end
        
        function list = getElement_soln(obj,list)
            if (isempty(obj.leftChild) && isempty(obj.rightChild))
%                 list=[list, obj.solution];
                list= [list, obj.solution(1)];
                
                
            else
                list = getElement_soln(obj.leftChild,list);
                list = getElement_soln(obj.rightChild,list);
            end
        end
        
        function list = getElement_bound(obj,list)
            if (isempty(obj.leftChild) && isempty(obj.rightChild))
%                 list{end+1}= obj.bound;
                list = [list, obj.bound(1)];
                
                
            else
                list = getElement_bound(obj.leftChild,list);
                list = getElement_bound(obj.rightChild,list);
            end
        end
        
        function list = getElements(obj, list)
            if (isempty(obj.leftChild) && isempty(obj.rightChild))
                list = [list, obj];
            else
                list = getElements(obj.leftChild,list);
                list = getElements(obj.rightChild,list);
            end
        end

        function list_outer = getElements_arr(obj, list)
            list_outer = {};
            
            if (isempty(obj.leftChild) && isempty(obj.rightChild))
                list = [list, obj];
            else
                list = getElements(obj.leftChild,list);
                list = getElements(obj.rightChild,list);
            end
            
        end
        function obj = assign_vals(obj, solution, bound)
            
            if (isempty(obj.leftChild) && isempty(obj.rightChild))
                obj.bound = bound;
                obj.solution= solution;
            else
                obj.leftChild = obj.leftChild.assign_vals(solution,bound);
                obj.rightChild = obj.rightChild.assign_vals(solution,bound);
                
            end
               
        end
        
    end
end