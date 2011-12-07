classdef Ellipsoid < handle
% Represents an E. coli bacteria in 3D space based on a collection of 2D
% ellipses.

    properties
        Center % [x, y, z]
        Orientation % [theta, phi]
        Size % [a, b] (assume both minor axes are the same size)
        Ellipses % vector of ellipses
        
        Claims = [];
    end
    
    properties (Dependent)
        theta, phi
        x, y, z
        maxZ
        a, b
        aspectRatio
    end
    
    properties (Constant)
        % TODO: rather arbitrary, but they seem to work well
        TYPICAL_ASPECT_RATIO = 6.0;
        LOW_ASPECT_RATIO = 1.20;
        MIN_SCORE = 0.5;
    end
    
    methods
        
        function obj = Ellipsoid(Ellipses)
            if nargin == 0, return, end
            
            obj.Ellipses = Ellipses;
            obj.recalculate();
        end
        
        function theta = get.theta(obj), theta = obj.Orientation(1); end
        function phi = get.phi(obj), phi = obj.Orientation(2); end
        function a = get.a(obj), a = obj.Size(1); end
        function b = get.b(obj), b = obj.Size(2); end
        function x = get.x(obj), x = obj.Center(1); end
        function y = get.y(obj), y = obj.Center(2); end
        function z = get.z(obj), z = obj.Center(3); end
        function z = get.maxZ(obj), z = max([obj.Ellipses.z]); end
        function r = get.aspectRatio(obj), r = obj.a/obj.b; end
        
        function recalculate(obj)
        % Recalculates the properties of the ellipsoid based on the current
        % Ellipses vector. The center is a weighted average of ellipse centers,
        % b is a weighted average of ellipse bs. a, theta, and phi each have
        % their own estimation functions.
        
            Areas = [obj.Ellipses.area];
            
            obj.Center = [weightedMean([obj.Ellipses.x], Areas), ...
                          weightedMean([obj.Ellipses.y], Areas), ...
                          weightedMean([obj.Ellipses.z], Areas)];
            
            obj.Size = [obj.estimateA(), weightedMean([obj.Ellipses.b], Areas)];
            
            obj.Orientation = [obj.estimateTheta(), obj.estimatePhi()];
        end
        
        function append(obj, ellipse)
        % Add an ellipse and recalculate.
        
            obj.Ellipses(end + 1) = ellipse;
            obj.recalculate();
        end
        
        function [Ellipses certainty] = calculateSearchAreas(obj, z)
        % Calculates a vector of ellipse-shaped regions to search for a possible
        % portion of this ellipsoid in zth x-y slice. Certainty is a measure of
        % how sure we are that the next ellipse will be there (higher is more
        % certain).
        
            certainty = length(obj.Ellipses);
            
            angle = obj.phi;
            Size = obj.Size;
            
            % If we've only found 1 ellipse so far, theta is ambiguous and so we
            % have 2 separate regions to search.
            if length(obj.Ellipses) == 1
                theta1 = obj.theta;
                theta2 = pi - theta1;
                dz = z - obj.z;
                
                r1 = dz/tan(theta1); r2 = dz/tan(theta2);
                C1 = [r1*cos(angle) + obj.x, r1*cos(angle) + obj.y, z];
                C2 = [r2*cos(angle) + obj.x, r2*cos(angle) + obj.y, z];
                
                Ellipses = [Ellipse(C1, Size, angle), Ellipse(C2, Size, angle)];
            else
            % Otherwise we can just extrapolate based on existing ellipses.
                X = [obj.Ellipses.x]; Y = [obj.Ellipses.y]; 
                Z = [obj.Ellipses.z];
                x = interp1(Z, X, z, 'linear', 'extrap');
                y = interp1(Z, Y, z, 'linear', 'extrap');
                
                Ellipses = Ellipse([x y z], Size, angle);
            end
        end
        
        function stakePossibleClaims(obj, Claims, matchEllipse)
        % Given the ellipse that was used to match claims, the ellipsoid will
        % place stakes in the given claims based on how well each claim matches
        % matchEllipse.
        
            Scores = zeros(size(Claims));
            for i = 1 : length(Claims)
                Scores(i) = ...
                    Ellipsoid.assignScore(Claims(i).object, matchEllipse);
            end
            
            Keep = Scores > obj.MIN_SCORE;
            Scores = Scores(Keep);
            Claims = Claims(Keep);
            
            % Rank preferences
            [Scores, Indices] = sort(Scores);
            
            numPrev = length(obj.Claims);
            obj.Claims = [obj.Claims Claims(Indices)];
            
            for i = 1 : length(Claims)
                Claims(i).stake(obj, Scores(i), i + numPrev);
            end
        end
        
        function revokeClaims(obj)
        % Revoke and clear all claims.
        
            for i = 1 : length(obj.Claims)
                obj.Claims(i).revoke(obj);
            end
            obj.Claims = [];
        end
        
    end
    
    methods (Hidden)
        
        function a = estimateA(obj)
        % Estimates teh length of the major axis (a) by finding the furthest
        % distance between two focal points from (possibly) different ellipses.
        % TODO: this isn't the best approximation, couldl be improved by
        % actually doing the math here, but for now its ok.
        
            if length(obj.Ellipses) == 1
                a = obj.Ellipses(1).a;
                return
            end
            
            furthest = 0;
            distSquared = @(P1, P2) sum((P1 - P2).^2);
            
            for i = 1:(length(obj.Ellipses) - 1)
                e1 = obj.Ellipses(i);
                for j = (i + 1):length(obj.Ellipses)
                    e2 = obj.Ellipses(j);
                    
                    furthest = max([furthest, distSquared(e1.F1, e2.F1), ...
                                              distSquared(e1.F1, e2.F2), ...
                                              distSquared(e1.F2, e2.F1), ...
                                              distSquared(e1.F2, e2.F2)]);
                end 
            end
            
            a = sqrt(furthest)/2;
        end
        
        function theta = estimateTheta(obj)
        % Tries to estimate theta -- the method changes based on how many
        % ellipses we've already seen.
        
            if length(obj.Ellipses) == 1
                theta = obj.estimateThetaUsingAspectRatio();
            else
                theta = obj.estimateThetaUsingPositions();
            end
        end
        
        function theta = estimateThetaUsingAspectRatio(obj)
        % Tries to measure theta by comparing the current aspect ratio against
        % the expected aspect ratio. This is relatively inaccurate, but its the
        % best we can do with 1 ellipse. NOTE: This is also ambiguous, the real
        % theta could be +/- the return result.
        
            ratio = mean([obj.Ellipses.aspectRatio]);
            theta = interp1([1, obj.TYPICAL_ASPECT_RATIO], [0, pi/2], ratio, ...
                            'linear', pi/2);
        end
        
        function theta = estimateThetaUsingPositions(obj)
        % Estimate theta by fitting a line to the x,z centers of the ellipses
        % and then measuring the slope.
        
            X = [obj.Ellipses.x]; X = X - X(1);
            if sum(X) < 1
                theta = 0.0;
                return
            end
            
            Z = [obj.Ellipses.z];
            P = polyfit(X - X(1), Z - Z(1), 1);
            theta = normalizeAngle(atan(1/P(1)));
        end
        
        function phi = estimatePhi(obj)
        % If there is more than one ellipse found, but they are roughly
        % circular, we'll use positions to estimate phi, otherwise we'll use the
        % angles of the ellipses in the x-y plane.
        
            meanRatio = mean([obj.Ellipses.aspectRatio]);
            if length(obj.Ellipses) > 1 && meanRatio < obj.LOW_ASPECT_RATIO
                phi = obj.estimatePhiUsingPositions();
            else
                phi = obj.estimatePhiUsingAngles();
            end
        end
        
        function phi = estimatePhiUsingPositions(obj)
        % Estimate phi by fitting a line to the x,y centers of the ellipses
        % and then measuring the slope.
        
            X = [obj.Ellipses.x]; X = X - X(1);
            if sum(X) < 1
                phi = pi/2;
                return
            end
            
            Y = [obj.Ellipses.y];
            P = polyfit(X, Y - Y(1), 1);
            phi = normalizeAngle(atan(P(1)));
        end
        
        function phi = estimatePhiUsingAngles(obj)
        % Take an average of the ellipse angles to estimate phi (ooh
        % non-euclidian statistics).
        
            Angles = [];
            for i = 1:length(obj.Ellipses)
                Angles(i) = normalizeAngle(obj.Ellipses(i).angle);
            end
            
            meanX = mean(cos(Angles));
            meanY = mean(sin(Angles));
            phi = normalizeAngle(atan2(meanY, meanX));
        end
        
    end
    
    methods (Static)
        function score = assignScore(ellipse, matchEllipse)
        % A simple scoring algorithim that assigns higher scores to better
        % matches based on the focal and size errors. TODO: improve this?
        
            score = 3/(1 + ellipse.getFocalError(matchEllipse)) + ...
                    1/(1 + ellipse.getSizeError(matchEllipse));
        end
    end
    
end

function v = weightedMean(Values, Weights)
    v = sum(Values .* Weights)/sum(Weights);
end

function angle = normalizeAngle(angle) 
% Between 0 and pi
    while angle > pi, angle = angle - pi; end
    while angle < 0, angle = angle + pi; end
end