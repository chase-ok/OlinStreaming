function f = velocityField(n, m)
    if nargin < 1, n = 10; end
    if nargin < 2, m = n; end
    
    
    function Result = calculate(imageSeries, TXYAdXdY)
        BinSize = imageSeries.ImageSize./[n, m];
        % Calculate center positions and initialze Velocities
        
        Positions = cell(n, m);
        Velocities = cell(n, m);
        Counts = zeros(n, m);
        for i = 1:n
            for j = 1:m
                Positions{i, j} = [i - 1, j - 1].*BinSize + BinSize/2;
                Velocities{i, j} = [0, 0];
            end
        end

        % Sum all of the velocities in each bin, keep tracking of count so that we
        % can average it later.
        for i = 1:size(TXYAdXdY, 1)
            % find position index based on x, y
            B = ceil(TXYAdXdY(i, 2:3)./BinSize);
            
            Velocities{B(1), B(2)} = Velocities{B(1), B(2)} + TXYAdXdY(i, 5:6);
            Counts(B(1), B(2)) = Counts(B(1), B(2)) + 1;
        end

        % finish average
        for i = 1:n
            for j = 1:m
                if Counts(i, j) == 0, continue, end % avoid 0 divide
                Velocities{i, j} = Velocities{i, j}/Counts(i, j);
            end
        end

        % flatten into matrices
        Result = {cell2mat(Positions(:)), cell2mat(Velocities(:))};
    end
    
    f = @calculate;
end