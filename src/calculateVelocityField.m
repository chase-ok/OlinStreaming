function [Positions Velocities] = calculateVelocityField(imageSeries, M, n, ...
                                  TimeRange)
% [Positions Velocities] = calculateVelocityField(imageSeries, M, n, TimeRange)
% Calculates a grided velocity field from analyzed particles paths over a given
% time range.
%
% Arguments:
%     - imageSeries: ImageSeries instance returned by analyzeParticles
%     - M: particle path matrix returned by analyzeParticles
%     - n: the number of bins to use per dimension. For example, n = 10 uses a
%     10-by-10 grid of bins to calculate the velocity field.
%     - TimeRange: a vector with two elements: (1) the earliest time to include
%     in the field, (2) the latest time to include. For example, TimeRange = [5
%     5.5] uses all particle positions recorded between 5 and 5.5 seconds.
%
% Returns:
%     - Positions: an n*n-by-2 matrix where each row is the center position (in
%     microns) of a bin. First column is x, second column is y.
%     - Velocities: an n*n-by-2 matrix where each row is the corresponding
%     average velocity of the bin in the Positions matrix. First column is Vx, 
%     second column is Vy.
% 
% Usage:
%     Current folder is streaming,
%     [imageSeries M] = analyzeParticles('data\July292011\High density', ...
%                                        80, 0, 0, true);
%     [P V] = calculateVelocityField(imageSeries, M, 32);
%     plotVelocityField(P, V);
% 
    BinSize = imageSeries.ImageSize/n;

    % Calculate center positions and initialze Velocities
    Positions = cell(n, n);
    Velocities = cell(n, n);
    Counts = zeros(n, n);
    for i = 1:n
        for j = 1:n
            Positions{i, j} = [i - 1, j - 1].*BinSize + BinSize/2;
            Velocities{i, j} = [0, 0];
        end
    end
    
    % Select only the data that falls within the given time range
    M = M(TimeRange(1) <= M(:, 1) & M(:, 1) <= TimeRange(2), :);
    
    % Sum all of the velocities in each bin, keep tracking of count so that we
    % can average it later.
    for i = 1:size(M, 1)
        B = ceil(M(i, 2:3)./BinSize); % find position index based on x, y
        Velocities{B(1), B(2)} = Velocities{B(1), B(2)} + M(i, 5:6);
        Counts(B(1), B(2)) = Counts(B(1), B(2)) + 1;
    end
    
    % finish average
    for i = 1:n
        for j = 1:n
            if Counts(i, j) == 0, continue, end % avoid 0 divide
            Velocities{i, j} = Velocities{i, j}/Counts(i, j);
        end
    end
    
    % flatten into matrices
    Positions = cell2mat(Positions(:));
    Velocities = cell2mat(Velocities(:));
end