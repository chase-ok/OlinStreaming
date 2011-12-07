function Objects = mergeEcoliSlices(Slices)
% TODO: REMOVE DEPENDENCE ON THIS FUNCTION -- (I think it's no longer used...)
% This is old and has been subsumed by the ZStackManager.

    MIN_SHARED_AREA_RATIO = 0.75;
    
    Used = makeUsed(Slices);
    
    for i_ = 1:length(Slices)
        for j_ = 1:length(Slices{i_})
            if Used{i_}(j_), continue, end
            
            [SliceIndices EllipseIndices] = group(i_, j_);
            parseGroup(SliceIndices, EllipseIndices);
        end
    end
    
    Objects = [];
    
    function [SliceIndices EllipseIndices] = group(initSlice, initEllipse)
        function MatchIndices = findInSlice(matchSlice, baseSlice, BaseIndices)
            MatchIndices = [];
            MatchSlice = Slices{matchSlice};
            BaseSlice = Slices{baseSlice};
            
            for i = 1:length(BaseIndices)
                ellip = BaseSlice(i);
                
                for j = 1:length(MatchSlice)
                    if Used{matchSlice}(j) || any(MatchIndices == j)
                        continue
                    end
                    
                    sharedArea = ellip.areaOfIntersection(MatchSlice(j));
                    if sharedArea == 0, continue, end

                    ratio = sharedArea/(min([ellip.area MatchSlice(j).area]));
                    if ratio > MIN_SHARED_AREA_RATIO
                        MatchIndices(end + 1) = j;
                    end
                end
            end
        end
        
        function [SliceIndices Groups] = search(SliceRange)
            SliceIndices = [];
            Groups = {};
            
            baseSlice = initSlice;
            BaseIndices = initEllipse;
            
            for i = SliceRange
                SliceIndices(end + 1) = i;
                Groups{end + 1} = findInSlice(i, baseSlice, BaseIndices);
                
                baseSlice = i;
                BaseIndices = Groups{end};
                if isempty(BaseIndices), break, end
            end
        end
        
        [LowerIndices LowerGroups] = search(initSlice:-1:1);
        [UpperIndices UpperGroups] = search((initSlice + 1):length(Slices));
        
        SliceIndices = [LowerIndices(end:-1:1) UpperIndices];
        EllipseIndices = [LowerGroups(end:-1:1) UpperGroups];
        markUsed(SliceIndices, EllipseIndices);
    end

    function markUsed(SliceIndices, EllipseIndices)
        for i = 1:length(SliceIndices)
            for j = 1:length(EllipseIndices)
                Used{i}(j) = true;
            end
        end
    end

    function parseGroup(SliceIndices, EllipseIndices)
        disp('--- GROUP ---');
        disp([num2str(length(SliceIndices)) ' ' num2str(length(EllipseIndices))]);
    end
end

function Used = makeUsed(Slices)
    Used = {};
    for i = 1:length(Slices)
        Used{i} = false(size(Slices{i}));
    end
end


% function Objects = mergeEcoliSlices(Slices)
% 
%     SEARCH_FACTOR = 1.5; % Search within an ellipse expanded by this amount
%     
%     Objects = {};
%     numObjects = 0;
%     Slices = addUsedField(Slices);
%     
%     numSlices = length(Slices);
%     for i = 1 : numSlices
%         Slices{i}.Orientations = Slices{i}.Orientations*pi/180; % to radians
%     end
%     
%     % search bottom to top, merging as we go
%     for sliceIndex = 1 : numSlices
%         disp(['Slice ' num2str(sliceIndex) ' of ' num2str(numSlices)]);
%         Slice = Slices{sliceIndex};
%         
%         for objId = 1 : length(Slice.X)
%             if Slices{sliceIndex}.Used(objId), continue, end
%             
%             Slices{sliceIndex}.Used(objId) = true;
%             
%             X = Slice.X(objId); Y = Slice.Y(objId); Z = Slice.Z(objId);
%             A = Slice.A(objId); B = Slice.B(objId);
%             Phi = Slice.Orientations(objId);
%             
%             for merge = (sliceIndex + 1) : numSlices
%                 % TODO: Maybe give it a 2 slice memory?
%                 MSlice = Slices{merge};
%                 
%                 x0 = X(end); y0 = Y(end);
%                 a0 = A(end); b0 = B(end);
%                 phi0 = Phi(end);
%                 
%                 [I Distances] = findObjectsNearby(MSlice, x0, y0, ...
%                                                   SEARCH_FACTOR*a0, ...
%                                                   SEARCH_FACTOR*b0, phi0);
%                 % TODO: maybe try enlarging search area
%                 if isempty(I), break, end
%                 
%                 Error = getError(Distances, a0, b0, phi0, MSlice.A(I), ...
%                                  MSlice.B(I), MSlice.Orientations(I));
%                 
%                 Sorted = sortrows([I Error], 2);
%                 best = Sorted(1, 1);
%                 
%                 X = [X MSlice.X(best)]; Y = [Y MSlice.Y(best)]; 
%                 Z = [Z MSlice.Z(best)];
%                 A = [A MSlice.A(best)]; B = [B MSlice.B(best)];
%                 Phi = [Phi MSlice.Orientations(best)];
%                 
%                 Slices{merge}.Used(best) = true;
%             end
%             
%             if length(X) >= 3
%                 numObjects = numObjects + 1;
%                 Objects{numObjects} = parseObject(X, Y, Z, A, B, Phi);
%             end
%         end
%     end
% end
% 
% function Error = getError(Distances, a0, b0, phi0, A, B, Phi)
%     DistanceError = Distances.^2;
%     PhiError      = (abs(Phi - phi0)/pi).^2;
%     SizeError     = (abs(A - a0)/a0).^2 + (abs(B - b0)/b0).^2;
%     
%     Error = 6*DistanceError + 3*PhiError + SizeError;
% end
% 
% function Slices = addUsedField(Slices)
%     for j = 1 : length(Slices)
%         Slices{j}.Used = false(size(Slices{j}.X));
%     end
% end
% 
% function v = weightedMean(Values, Weights)
%     v = sum(Values .* Weights)/sum(Weights);
% end
% 
% function Object = parseObject(X, Y, Z, A, B, Phi)
%     ECCENTRICITY_THRESHOLD = 2.0;
%     
%     Areas = pi*A.*B;
%     Object.x = weightedMean(X, Areas);
%     Object.y = weightedMean(Y, Areas);
%     Object.z = weightedMean(Z, Areas);
%     Object.b = weightedMean(B, Areas);
%     
%     Object.numSlices = length(X);
%     
%     [a, aIndex] = max(A);
%     Object.a = a;
%     
%     if Object.numSlices <= 3
%         Object.phi   = weightedMean(Phi, Areas);
%         Object.theta = pi/2;
%     else
%         if Object.a/Object.b > ECCENTRICITY_THRESHOLD
%             Object.phi = Phi(aIndex);
%         else
%             P = polyfit(X - X(1), Y - Y(1), 1);
%             Object.phi = atan(P(1)); % TODO: this might need to be negative...
%         end
%         
%         P = polyfit(X - X(1), Z - Z(1), 1);
%         Object.theta = atan(1/P(1));
%     end
% end