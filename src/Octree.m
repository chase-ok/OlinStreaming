classdef Octree < handle
% Recursively splits a 3D region into 8 smaller regions.
% TODO: no longer necessary?

    properties
        TopLeft
        Size
        Objects
    end
    
    properties (Dependent)
        numObjects
    end
    
    methods
        
        function obj = Octree(TopLeft, Size, ObjectPool)
            if nargin == 0, return, end
            
            obj.TopLeft = TopLeft;
            obj.Size = Size;
            obj.select(ObjectPool);
        end
        
        function num = get.numObjects(obj)
            num = length(obj.Objects);
        end
        
        function select(obj, ObjectPool)
            BottomRight = obj.TopLeft + obj.Size;
            
            for i = 1:length(ObjectPool)
                test = ObjectPool(i);
                Inside(i) = all(test.Center > obj.TopLeft) && ... 
                            all(test.Center < BottomRight);
            end
            
            obj.Objects = ObjectPool(Inside);
        end
        
        function Trees = divide(obj)
            NewSize = obj.Size/2;
            make = @(offset) Octree(obj.TopLeft + offset, NewSize, obj.Objects);
            
            Trees(1) = make([0          0          0]);
            Trees(2) = make([NewSize(1) 0          0]);
            Trees(3) = make([0          NewSize(2) 0]);
            Trees(4) = make([NewSize(1) NewSize(2) 0]);
            Trees(5) = make([0          0          NewSize(3)]);
            Trees(6) = make([NewSize(1) 0          NewSize(3)]);
            Trees(7) = make([0          NewSize(2) NewSize(3)]);
            Trees(8) = make([NewSize(1) NewSize(2) NewSize(3)]);
        end
    end
    
end

