classdef ZStackManager < handle
% Handles merge ellipses found in x,y slices into 3D ellipsoids.

    properties
        Size % [x, y]
        RegionSize % [x, y]
        Regions % NUM_REGIONS-by-NUM_REGIONS cell matrix
        zStep % z-difference between each x,y slice
    end
    
    properties (Hidden)
        PendingEllipsoids % Ellipsoids that are still looking for matches
        Ellipsoids % Complete ellipsoids
        currentZ
        debug = false;
    end
    
    properties (Constant)
        NUM_REGIONS = [16 16];
        EXPANSION = 0.3; % How much to expand search areas by
        CERTAINTY_ADDITION = 2; % Certainty adjustments
        MEMORY = 2.0 + 0.0001; % Number of slices to "remember" an ellipsoid
    end
    
    methods
        function obj = ZStackManager(width, height, zStep, debug)
            if nargin == 0, return, end
            if nargin == 3, debug = false; end
            
            obj.Size = [width height];
            obj.RegionSize = obj.Size ./ obj.NUM_REGIONS;
            obj.zStep = zStep;
            obj.debug = debug;
        end
        
        function Ellipsoids = findEllipsoids(obj, Slices)
        % Main function -- Takes a set of x,y slices and returns all of the
        % ellipsoids found in them as a vector.
        
            obj.Ellipsoids = Ellipsoid.empty();
            obj.PendingEllipsoids = Ellipsoid.empty();
            
            for i = 1:length(Slices)
                obj.print(['Finding ellipsoids in slice #' num2str(i)]);
                obj.currentZ = (i - 1)*obj.zStep;
                obj.parseSlice(Slices{i});
            end
            
            Ellipsoids = [obj.Ellipsoids obj.PendingEllipsoids];
        end
        
    end
    
    methods (Hidden)
        
        function parseSlice(obj, Slice)
        % Called on each slice in order (accending z)
        
            obj.createClaimRegions(Slice);
            obj.print('Claim regions created.');
            obj.stakeClaims();
            obj.print('Claims staked');
            Unresolved = obj.resolveClaims();
            obj.print('Claims resolved');
            obj.purgePendingEllipsoids();
            obj.print('Pending purged');
            obj.createNewEllipsoids(Unresolved);
            obj.print('New ellipses created');
        end
        
        function createClaimRegions(obj, Slice)
        % Create a claim for each ellipse in this slice and then divide the
        % claims up into regions so that they can be found faster.
        
            obj.Regions = cell(obj.NUM_REGIONS);
            RS = obj.RegionSize; % faster
            
            for i = 1:length(Slice)
                ellip = Slice(i);
                XY = ceil(ellip.Center(1:2)./RS);
                obj.Regions{XY(2), XY(1)}(end + 1) = Claim(ellip);
            end
        end
        
        function stakeClaims(obj)
        % For each pending ellipsoid, calculate a search area for matching
        % ellipses, and then stake claims on them based on how good of a match
        % they are.
        
            n = length(obj.PendingEllipsoids);
            for i = 1:n
                ellip = obj.PendingEllipsoids(i);
                
                [Ellipses cert] = ellip.calculateSearchAreas(obj.currentZ);
                for j = 1:length(Ellipses)
                    Claims = obj.getClaims(Ellipses(j), cert);
                    ellip.stakePossibleClaims(Claims, Ellipses(j));
                end
                
                obj.print(['Claim staked (' num2str(i) '/' num2str(n) ')']);
            end
        end
        
        function Claims = getClaims(obj, ellipse, cert)
        % Given a degree of certainity and an ellipse-shaped search area, find
        % all possible matching ellipses in the current slice.
        
            % Expand the search area if we are less certain
            a = ellipse.a*(1 + sqrt(ellipse.a)* ...
                   (obj.EXPANSION + obj.CERTAINTY_ADDITION / (1 + cert)));
            aSquared = a^2;
            
            % Check the claim regions that fall within the reach of the search
            % area
            Offsets = { [0, 0], [a, a], [-a, a], [a, -a], [-a, -a] };
            Center = ellipse.Center(1:2);
            RS = obj.RegionSize;
            XYs = {};
            for i = 1:length(Offsets)
                Pos = Center + Offsets{i};
                XY = [min(obj.NUM_REGIONS(1), max(1, ceil(Pos(1)/RS(1)))), ...
                      min(obj.NUM_REGIONS(2), max(1, ceil(Pos(2)/RS(2))))];
                
                found = false;
                for j = 1:length(XYs)
                    if all(XY == XYs{j})
                        found = true;
                        break
                    end
                end
                if ~found, XYs{end + 1} = XY; end
            end
            
            % Filter out too-distant claims
            Claims = Claim.empty();
            for i = 1:length(XYs)
                Region = obj.Regions{XYs{i}(2), XYs{i}(1)};
                
                for j = 1:length(Region)
                    match = Region(j);
                    
                    distSquared = sum((Center - match.object.Center(1:2)).^2);
                    if distSquared < aSquared
                        Claims(end + 1) = match;
                    end
                end
            end
        end
        
        function Unresolved = resolveClaims(obj)
        % Go through all of the claims and pick a winning ellipsoid and assign
        % the claim to the winner. Any claims that don't have a winner become
        % the start of a new ellipsoid.
        
            Unresolved = Claim.empty();
            
            total = num2str(numel(obj.Regions));
            n = 1;
            for i = 1:size(obj.Regions, 1)
                for j = 1:size(obj.Regions, 2)
                    obj.print(['Resolving claims in region ' num2str(n), ...
                               ' of ' total]);
                    n = n + 1;
                    Region = obj.Regions{i, j};
                    
                    for k = 1:length(Region)
                        claim = Region(k);
                        
                        if claim.hasStakes
                            ellipsoid = claim.getWinner();
                            ellipsoid.revokeClaims();
                            ellipsoid.append(claim.object);
                        else
                            Unresolved(end + 1) = claim;
                        end
                    end
                end
            end
        end
        
        function purgePendingEllipsoids(obj)
        % If an ellipsoid hasn't found a match in a while, purge it from the
        % pending list and add it to the complete ellipsoid list.
        
            zLimit = obj.currentZ - obj.MEMORY*obj.zStep;
            Keep = [obj.PendingEllipsoids.maxZ] > zLimit;
            obj.Ellipsoids = [obj.Ellipsoids obj.PendingEllipsoids(~Keep)];
            obj.PendingEllipsoids = obj.PendingEllipsoids(Keep);
        end
        
        function createNewEllipsoids(obj, Claims)
            New(length(Claims)) = Ellipsoid; % preallocate
            for i = 1:length(Claims)
                New(i) = Ellipsoid(Claims(i).object);
            end
            obj.PendingEllipsoids = [obj.PendingEllipsoids New];
        end
        
        function print(obj, msg)
            if obj.debug
                disp(msg)
            end
        end
    end
    
end

