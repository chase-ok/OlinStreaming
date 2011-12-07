classdef Quadtree < handle
% Recursively divide a 2D region into 4 smaller regions.
% TODO: no longer needed?
    
    properties
        TopLeft
        Size
        Objects
    end
    
    properties (Constant)
        BUFFER = 0.05
    end
    
    properties (Dependent)
        numObjects
    end
    
    methods
        function obj = Quadtree(TopLeft, Size, ObjectPool)
            if nargin == 0, return, end
            
            obj.TopLeft = TopLeft;
            obj.Size = Size;
            obj.select(ObjectPool);
        end
        
        function num = get.numObjects(obj)
            num = length(obj.Objects);
        end
        
        function select(obj, ObjectPool)
            TopLeft = obj.TopLeft - obj.Size*obj.BUFFER;
            BottomRight = obj.TopLeft + obj.Size*(1 + obj.BUFFER);
            
            for i = 1:length(ObjectPool)
                test = ObjectPool(i);
                Inside(i) = all(test.Center(1:2) > TopLeft) && ... 
                            all(test.Center(1:2) < BottomRight);
            end
            
            obj.Objects = ObjectPool(Inside);
        end
        
        function Trees = subdivide(obj)
            NewSize = obj.Size/2;
            make = @(p) Quadtree(obj.TopLeft + p, NewSize, obj.Objects);
            
            Trees(1) = make([0          0         ]);
            Trees(2) = make([NewSize(1) 0         ]);
            Trees(3) = make([0          NewSize(2)]);
            Trees(4) = make([NewSize(1) NewSize(2)]);
        end
    end
    
    methods (Static)
        function Trees = divide(Objects, width, height, n)
            Trees = Quadtree([0 0], [width height], Objects);
            for i = 1 : n
                NewTrees = [];
                for j = 1 : length(Trees)
                    NewTrees = [NewTrees Trees(i).subdivide()];
                end
                Trees = NewTrees;
            end
        end
        
        function TreeSlices = convertSlices(Slices, width, height, n)
            for i = 1 : length(Slices)
                TreeSlices{i} = Quadtree.divide(Slices{i}, width, height, n);
            end
        end
    end
    
end

