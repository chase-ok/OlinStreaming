function calculate = byRadius(f, shouldAverage, dr, maxRadius)
    if nargin == 2
        Radiuses = [0.5:0.25:8, 9:20];
    elseif nargin == 3
        Radiuses = dr;
    elseif nargin == 4
        Radiuses = 0.5:dr:maxRadius;
    else
        error('Too many arguments.');
    end

    RadiusesSq = Radiuses.^2;
    Areas = pi.*[RadiusesSq(1), (RadiusesSq(2:end) - RadiusesSq(1:end-1))];
    
    function Result = wrapper(imageSeries, TXYAdXdY)
        XY = TXYAdXdY(:, 2:3);
        n = size(XY, 1);
        
        Data = zeros(size(Radiuses));
        
        for i = 1:n
            DistancesSq = sum((XY - repmat(XY(i, :), n, 1)).^2, 2);
            
            for j = 1:length(RadiusesSq)
                if j == 1
                    rSq1 = 0.001; % Don't include the measurment particle
                else
                    rSq1 = RadiusesSq(j-1);
                end
                rSq2 = RadiusesSq(j);
                Rows = (DistancesSq >= rSq1) & (DistancesSq <= rSq2);
                
                Data(j) = Data(j) + f(TXYAdXdY(i, :), TXYAdXdY(Rows, :));
            end
        end
        
        if shouldAverage
            n = sum(DistancesSq <= RadiusesSq(end));
            Data = Data./n;
        end
        Density = Data./Areas;
        Result = {Radiuses, Density};
    end
    
    calculate = @wrapper;
end