classdef Ellipse
% Represents an ellipse shape detected in an image.

    properties
        Center % [x, y, z] (in some contexts z might be t)
        Size   % [a, b] (long and short axes, respectively)
        angle  % radians 
        
        used = false; % used in a 3D bacteria group
    end
    
    properties (Dependent)
        x, y, z
        a, b
        area
        
        eccentricity
        aspectRatio % Ratio of a to b
        
        focalDistance
        F1, F2 % focal points [x, y]
    end
    
    methods
        function obj = Ellipse(Center, Size, angle)
            if nargin == 0, return, end
            
            obj.Center = Center;
            obj.Size   = Size;
            obj.angle  = angle;
        end
        
        function a = get.a(obj), a = obj.Size(1); end
        function b = get.b(obj), b = obj.Size(2); end
        function x = get.x(obj), x = obj.Center(1); end
        function y = get.y(obj), y = obj.Center(2); end
        function z = get.z(obj), z = obj.Center(3); end
        
        function area = get.area(obj)
            area = pi*prod(obj.Size);
        end
        
        function ecc = get.eccentricity(obj)
            ecc = sqrt(1 - (obj.b/obj.a)^2);
        end
        
        function f = get.focalDistance(obj)
            f = sqrt(obj.a^2 - obj.b^2);
        end
        
        function F1 = get.F1(obj)
            F1 = obj.Center(1:2) ...
                 + obj.focalDistance*[cos(obj.angle) sin(obj.angle)];
        end
        
        function F2 = get.F2(obj)
            F2 = obj.Center(1:2) ...
                 - obj.focalDistance*[cos(obj.angle) sin(obj.angle)];
        end
        
        function r = get.aspectRatio(obj)
            r = obj.a/obj.b;
        end
        
        function area = areaOfIntersection(obj, other)
        % Returns the area of the overlapping regions of two ellipses. NOTE:
        % ignores the z-component of the ellipse positions.
            
            % make relative to obj.Center
            other.Center = other.Center - obj.Center;
            obj.Center = zeros(size(obj.Center));
            
            distSquared = sum(other.Center.^2);
            if distSquared > (obj.a + other.a)^2
                area = 0;
                return
            end
            
            % now shift both ellipses so that everything is in the first
            % quadrant
            OffsetXY = -min([obj.Center(1:2) - obj.a 
                             other.Center(1:2) - other.a]);
            obj.Center(1:2)   = OffsetXY;
            other.Center(1:2) = other.Center(1:2) + OffsetXY;
            
            % make a blank, binary image to draw on
            Canvas = false(max([obj.Center(1:2) + obj.a
                                other.Center(1:2) + other.a]));
            Intersection = obj.fill(Canvas) & other.fill(Canvas);
            area = bwarea(Intersection);
        end
        
        function Image = fill(obj, Image)
        % Draws a filled-in ellipse to a binary image, overriding any pixels
        % underneath it.
        
            twoA = 2*obj.a;
            ImageSize = size(Image);
            
            minX = round(max([1, obj.x - obj.a])); maxX = round(min([ImageSize(2), obj.x + obj.a]));
            minY = round(max([1, obj.y - obj.a])); maxY = round(min([ImageSize(1), obj.y + obj.a]));
            
            % mesh of all possible pixels touched by the ellipse
            [X Y] = meshgrid(minX:maxX, minY:maxY);
            X = X(:); Y = Y(:); % flatten
            
            % use distance formula to determine which pixels are actually
            % touched by the ellipse
            Inside = twoA > sqrt((X - obj.F1(1)).^2 + (Y - obj.F1(2)).^2) + ...
                            sqrt((X - obj.F2(1)).^2 + (Y - obj.F2(2)).^2);
            X = X(Inside);
            Y = Y(Inside);
            
            % fill in the pixels
            for i = 1:length(X)
                Image(X(i), Y(i)) = true;
            end
        end
        
        function plot(obj, style, includeX)
        % Plots the ellipse shape using the given line style. If includeX is
        % unspecified or true, it also draws a red 'x' at the center of the
        % ellipse.
        
            angle = -obj.angle;
            T = linspace(0, 2*pi, 36);
            
            X = obj.x + (obj.a*cos(T)*cos(angle) - obj.b*sin(T)*sin(angle));
            Y = obj.y + (obj.a*cos(T)*sin(angle) + obj.b*sin(T)*cos(angle));
            plot(X, Y, style, 'LineWidth', 2);
            
            if nargin < 3 || includeX
                plot(obj.x, obj.y, 'rx');
            end
        end
        
        function Ellipses = split(obj, newArea, axis)
        % Splits an ellipse into multiple new ellipses each with (approximately)
        % newArea. Axis (either 'major' or 'minor') is the axis to split the
        % ellipse on. Used by splitLargeEllipses.
        
            switch axis
                case 'major'
                    a = obj.a;
                    b = obj.b;
                    angle = obj.angle;
                    
                case 'minor'
                    a = obj.b;
                    b = obj.a;
                    angle = obj.angle + pi/2;
            end
            
            n = round(obj.area/newArea);
            Ellipses(n) = Ellipse; % pre-allocate
            
            spacing = 2*a/(n + 1);
            offset = -a;
            
            for i = 1 : n
                Center = [obj.x + offset*cos(angle), ...
                          obj.y + offset*sin(angle), obj.z];
                Size = [b, spacing/2];
                Ellipses(i) = Ellipse(Center, Size, angle);
            end
        end
        
        function err = getFocalError(obj, other)
        % Returns a measurement of the focal point error between this ellipse 
        % and a possible matching ellipse. Essentially it returns the sum of the
        % distances between the matching focal points of each ellipse, the
        % greater this distance, the less 'similar' these two ellipses are.
        
            distSquared = @(P1, P2) sum((P1 - P2).^2);
            
            f1f1 = distSquared(obj.F1, other.F1);
            f1f2 = distSquared(obj.F1, other.F2);
            f2f1 = distSquared(obj.F2, other.F1);
            f2f2 = distSquared(obj.F2, other.F2);
            
            err = min([f1f1 + f2f2, f1f2 + f2f1]);
        end
        
        function err = getSizeError(obj, other)
        % Returns a measurement of the size error between this ellipse and a
        % possible matching ellipse. Returns 0 if they have the same long-axis.
        % The greater the return value, the more different they are in size.
        
            if obj.a > other.a
                a = obj.a;
                b = other.a;
            else
                a = other.a;
                b = obj.a;
            end
            
            err = (a - b)/a;
        end
    end
    
    methods (Static)
        function Ellipses = fromRegionProps(Props, z)
        % Takes the results of a call to regionprops and builds a vector of
        % ellipses from it. regionprops must have been called with the
        % appropriate arguments so that the Centroid, MajorAxisLength,
        % MinorAxisLength, and Orientaton properties are measured.
        
            n = length(Props);
            Ellipses(n) = Ellipse; % pre-allocate
            
            for i = 1 : n
                p = Props(i);
                
                Center = [p.Centroid, z];
                Size   = [p.MajorAxisLength/2, p.MinorAxisLength/2];
                angle  = p.Orientation * (pi/180);
                Ellipses(i) = Ellipse(Center, Size, angle);
            end
        end
        
        function Ellipses = splitLargeEllipses(Ellipses)
        % Analyzes the size distribution of a given set of ellipses, picks out
        % the largest outliers and splits them into smaller ellipses. Usually
        % this works well to break apart intersecting/touching ellipses in the
        % image, but sometimes runs into problems with confocal artifacts.
        
            ECCENTRICITY_THRESHOLD = 0.2;
            MAX_SIZE_FACTOR        = 2.0;
            OUTLIER_FACTOR         = 4.0;

            % Sort by eccentricities to determine which axis to split on
            % Also, size distribution is heavily dependent on eccentricity so
            % we'll split it up into 2 separate distributions
            
            Eccs = [Ellipses.eccentricity];
            LowEccs = Ellipses(Eccs < ECCENTRICITY_THRESHOLD);
            HighEccs = Ellipses(Eccs >= ECCENTRICITY_THRESHOLD);
            
            % Find outliers based on median measurements
            LowMedian = median([LowEccs.area]) * OUTLIER_FACTOR;
            HighMedian = median([HighEccs.area]) * OUTLIER_FACTOR;
            
            LowSplits = [LowEccs.area] > LowMedian*MAX_SIZE_FACTOR;
            HighSplits = [HighEccs.area] > HighMedian*MAX_SIZE_FACTOR;
            
            % Keep ellipses that don't need to be split
            OldEllipses = [LowEccs(~LowSplits) HighEccs(~HighSplits)];
            NewEllipses = [];
            
            ToSplit = LowEccs(LowSplits);
            for i = 1:length(ToSplit)
                ellip = ToSplit(i);
                NewEllipses = [NewEllipses ellip.split(LowMedian, 'major')];
            end
            
            ToSplit = HighEccs(HighSplits);
            for i = 1:length(ToSplit)
                ellip = ToSplit(i);
                NewEllipses = [NewEllipses ellip.split(HighMedian, 'minor')];
            end
            
            Ellipses = [OldEllipses NewEllipses];
        end
        
        function Image = drawBW(Ellipses, width, height)
        % Draws a vector of ellipses to a width-by-height binary image. 
            Image = false(height, width);
            for i = 1:length(Ellipses)
                Image = Ellipses(i).fill(Image);
            end
        end
        
        function M = toTrackingMatrix(Ellipses)
        % Turns an ellipse vector into a matrix that can be read by the track.m
        % algorithm. M has the following columns: (1) x, (2) y, (3) id, (4),
        % t/z. id is just the index of the ellipse inside the given ellipse
        % vector and is included so that we can find it again later.
        
            M = zeros(length(Ellipses), 4);
            M(:, 1) = [Ellipses.x];
            M(:, 2) = [Ellipses.y];
            M(:, 3) = 1:length(Ellipses); % id
            M(:, 4) = [Ellipses.z];
        end
    end
    
end

