function Frames = calculateFrameByFrame(func, imageSeries, TXYAdXdY, ...
                                        TimeRange, window)
    if nargin < 5, window = 1; end
    if nargin < 4, TimeRange = []; end
    
    if mod(window, 2) ~= 1
        error('Window must be a positive, odd integer (defaults to 1).');
    end
    
    Times = TXYAdXdY(:, 1);
    
    if ~isempty(TimeRange)
        TXYAdXdY = TXYAdXdY(TimeRange(1) <= Times & Times <= TimeRange(2), :);
    end
    
    % Break into frames

    FrameRows = zeros(0, 2);
    startRow = 1;
    frameTime = Times(1);
    for row = 1:size(TXYAdXdY, 1)
        rowTime = Times(row);
        
        if abs(frameTime - rowTime) > 0.00001
            FrameRows(end + 1, :) = [startRow, row - 1];
            
            frameTime = rowTime;
            startRow = row;
        end
    end
    
    disp(length(FrameRows));
    
    Frames = {};
    offset = (window - 1)/2;
    
    for i = (1 + offset):(length(FrameRows) - offset)
        startRow = FrameRows(i - offset, 1);
        endRow = FrameRows(i + offset, 2);
        FrameData = TXYAdXdY(startRow:endRow, :);
        Frames{end + 1} = func(imageSeries, FrameData);
    end
    
end