function [imageSeries TXYAdXdY] = analyzeParticles(folder, seriesNum, ...
                                  channelNum, timeOffset, debug)
% [imageSeries TXYAdXdY] = analyzeParticles(folder, series, channel, ...
%                                           timeOffset, debug)
% Calculates and analyzes the particle paths in a given image series. Assumes
% tracking E. coli bacteria (for now).
%
% Arguments:
%     - folder: string containing the relative path to the folder containing the
%     image series.
%     - seriesNum: integer series number
%     - channelNum: integer channel number
%     - timeOffset (defaults to 0): the starting time used for the series
%     (seconds). For example, if the series was taken 3 minutes after the sample
%     was prepared, pass 3*60 as the time offset.
%     - debug (defaults to false): if true, extra debugging information is
%     printed to the screen.
%
% Returns:
%     - imageSeries: ImageSeries instance describing the setting used to image
%     the sample, the number of frames, starting time, etc.
%     - TXYAdXdY: matrix with 1 row for each particle position and the following
%     collumns (in order): time (sec), x (microns), y (microns), angle
%     (radians), dx (x-velocity, microns/s), and dy (y-velocity, microns/s).
% 
% Usage:
%     Current path is set to the streaming folder,
%     [imageSeries, M] = analyzeParticles('data\July292011\High density', ...
%                                         80, 0, 0, true);
% 

    if nargin < 5, debug = false; end
    if nargin < 4, timeOffset = 0; end
    
    imageSeries = ImageSeries(folder, seriesNum, channelNum);
    imageSeries.load();
    Paths = trackEcoli(imageSeries.Images, debug);
    imageSeries.releaseImages();
    
    % tracking algorithm requires integer time intervals and matlab uses pixel
    % measurements for size / position so we have to correct all of the units
    % afterwards.
    for i = 1:length(Paths)
        Paths(i).correctUnits(timeOffset, imageSeries.dt, imageSeries.Pixel);
    end
    
    M = [];
    for path = Paths
        M = [M; path.getData()];
    end
    TXYAdXdY = sortrows(M, 1); % sort by time
end