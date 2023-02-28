classdef NewCylinder
    properties(Dependent)
        R
        Height
    end
    properties(Access=private)
        R_
        Height_
    end
    methods
        function C = NewCylinder(R, Height)
            if nargin > 0
                if ~isa(R, 'double') || ~isa(Height, 'double')
                    error('R and Height must be double.');
                end
                
                if ~isequal(size(R), size(Height))
                    error('Dimensions of R and Height must match.');
                end
                for k = numel(R):-1:1
                    C(k).R_ = R(k);
                    C(k).Height_ = Height(k);
                end
            end            
        end
        
        function V = volume(C)
            V = pi .* [C.R_].^2 .* [C.Height_];
        end
        
        function C = set.R(C, R)
            checkValue(R);
            C.R_ = R;
        end

        function R = get.R(C)
            R = C.R_;
        end
        
        function C = set.Height(C, Height)
            checkValue(Height);
            C.Height_ = Height;
        end
        
        function Height = get.Height(C)
            Height = C.Height_;
        end        
    end
end

function checkValue(x)
    if ~isa(x, 'double') || ~isscalar(x)
        error('value must be a scalar double.');
    end
end