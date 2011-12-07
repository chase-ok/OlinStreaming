function Objects = mergeCyanoSlices(Slices)
% TODO: REMOVE DEPENDENCE ON THIS FUNCTION -- (I think it's no longer used...)
% This is old and has been subsumed by the ZStackManager, since spheres are just
% a special case of ellipsoids.

    Objects = cell();
    numObjects = 0;
    Slices = addUsedField(Slices);
    
    % search bottom to top, merging as we go
    numSlices = length(Slices);
    for sliceIndex = 1 : numSlices
        Slice = Slices{sliceIndex};
        
        for objId = 1 : length(Slice.X)
            if Slice.Used(objId)
                continue
            end
            Slice.Used(objId) = true;
            
            X = Slice.X(objId);
            Y = Slice.Y(objId);
            Z = Slice.Z(objId);
            Areas = Slice.Area(objId);
            
            for mergeIndex = (objId + 1) : numSlices
                found = false; % TODO: Maybe give it a 2 slice memory?
                
                % TODO: find all matches within a certain centroid tolerance,
                % then find the closest one.
                
                for matchId = 1 : length(Slices{mergeIndex}.X)
                    if Slices{mergeIndex}.Used(matchId)
                        continue
                    end
                    
                    matchX = Slices{mergeIndex}.X(matchId);
                    matchY = Slices{mergeIndex}.Y(matchId);
                    dx = matchX - X(end);
                    dy = matchY - Y(end);
                    
                    if dx^2 + dy^2 < Areas(end)/pi
                        X = [X matchX];
                        Y = [Y matchY];
                        Z = [Z Slices{mergeIndex}.Z(matchId)];
                        Areas = [Areas Slices{mergeIndex}.Areas(matchId)];
                        
                        Slices{mergeIndex}.Used = true;
                        found = true;
                        break
                    end
                end
                
                if ~found
                    break
                end
            end
            
            numObjects = numObjects + 1;
            Objects{numObjects} = parseObject(X, Y, Z, Areas);
        end
    end
end

function Slices = addUsedField(Slices)
    for j = 1 : length(Slices)
        Slices{j}.Used = false(size(Slices{j}.X));
    end
end

function v = weightedMean(Values, Weights)
    v = sum(Values .* Weights)/sum(Weights);
end

function Object = parseObject(X, Y, Z, Areas)
    Object.x = weightedMean(X, Areas);
    Object.y = weightedMean(Y, Areas);
    Object.z = weightedMean(Z, Areas);
    Object.r = sqrt(max(Areas)/pi); % TODO: maybe take z-range into account?
end