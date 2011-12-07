classdef EllipsePath < handle
% Represents the path of an ellipse-shaped particle over time.

    properties
        Ellipses % n vector of recorded ellipses
        Velocities % (n - 1)-by-2 matrix of 
        Times % n vector of time stamps that the ellipses were recorded on
    end
    
    properties (Dependent)
        MeanCenter % Average [x, y] position
        MeanVelocity % Average [x, y] velocity
        num % number of recorded ellipses
        minTime % earliest timestamp
        maxTime % latest timestamp
    end
    
    methods
        function obj = EllipsePath(Ellipses, Times)
            if nargin == 0, return; end
            
            obj.Ellipses = Ellipses;
            obj.Times = Times;
            obj.calculateVelocities();
        end
        
        function num = get.num(obj)
            num = length(obj.Ellipses);
        end
        
        function m = get.MeanCenter(obj)
            m = mean(vertcat(obj.Ellipses.Center));
            m = m(1:2);
        end
        
        function m = get.MeanVelocity(obj)
            m = mean(obj.Velocities);
        end
        
        function t = get.minTime(obj)
            t = min(obj.Times);
        end
        
        function t = get.maxTime(obj)
            t = max(obj.Times);
        end
        
        function ok = hasTime(obj, t)
        % Returns true iff the given time falls within (minTime, maxTime)
        
            ok = obj.minTime <= t && t <= obj.maxTime;
        end
        
        function V = getVelocity(obj, t)
        % Returns the approximate velocity at a given time t. For now, there is
        % a simple linear interpolation for t values between time stamps.
        
            V = interp1(obj.Times(1:(end - 1)), obj.Velocities, t, ...
                        'linear', 'extrap');
        end
        
        function e = getEllipse(obj, t)
        % Returns the ellipse recorded closest to time t.
        
            index = find(obj.Times <= t, 1, 'last');
            e = obj.Ellipses(index);
        end
        
        function D = getData(obj, t)
        % Returns a matrix with the following columns: (1) t, (2) x, (3) y, (4)
        % angle, (5) Vx, (6) Vy. If t is specified, the matrix only has 1 row
        % and contains data from the ellipses recorded closest to that time.
        % Otherwise, there are multiple rows containing all of the recorded
        % ellipse positions.
        
            if nargin < 2
                for i = 1:obj.num
                    e = obj.Ellipses(i);
                    if i == obj.num
                        V = obj.Velocities(i - 1, :);
                    else
                        V = obj.Velocities(i, :);
                    end
                    D(i, :) = [obj.Times(i), e.x, e.y, e.angle, V];
                end
            else
                e = obj.getEllipse(t);
                D = [t, e.x, e.y, e.angle, obj.getVelocity(t)];
            end
        end
        
        function correctUnits(obj, timeOffset, dt, Pixel)
        % Uses an initial time offset and time and space units to convert from
        % frame number and pixel units to real world numbers. If the pixels are
        % not equally scaled, the size and angles of the ellipses could be 
        % distorted.
        
            for i = 1:obj.num
                e = obj.Ellipses(i);
                
                t = timeOffset + obj.Times(i)*dt;
                obj.Times(i) = t;
                
                Center = [e.Center(1:2).*Pixel, t];
                Size = e.Size*Pixel(1); % Assuming square field of view
                obj.Ellipses(i) = Ellipse(Center, Size, e.angle);
            end
            
            obj.calculateVelocities();
        end
    end
    
    methods (Hidden)
        function calculateVelocities(obj)
        % Calculates the velocities of the ellipses using a simple derivative
        % approximation. TODO: use a more accurate approximation?
        
            obj.Velocities = [];
            for i = 1:(obj.num - 1)
                e1 = obj.Ellipses(i);
                e2 = obj.Ellipses(i + 1);
                
                obj.Velocities(i, :) = (e2.Center(1:2) - e1.Center(1:2))/...
                                       (obj.Times(i + 1) - obj.Times(i));
            end
        end
    end
end