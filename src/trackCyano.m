function Paths = trackCyano(Images, debug)
% Paths = trackCyano(Images)
% Tracks cyanobacteria over a t-series of images using the tracking algorithm in
% track.m. The particles paths are grouped into instances of EllipsePath for
% further analysis.
%
% Arguments:
%     - Images: a cell vector of Images.
%     - debug: a boolean that, when true, causes extra debugging information to
%     printed.
%
% Returns:
%     - Paths: A vector of Path instances that describe the paths of individual
%     bacteria.
% 
% Usage:
%     Paths = trackCyano(Images, true);

    function print(msg), if debug, disp(msg), end, end

    print('Finding bacteria');
    TimeSlices = {}; % Cell vector of particle positions, grouped by time
    for i = 1:length(Images)
        print([num2str(round(100*i/length(Images))) '%']);
        TimeSlices{i} = findCyanoBacteria(Images{i}, i, false);
    end

    % Format particle positions so that 
    TrackingMatrix = [];
    for i = 1:length(TimeSlices)
        TrackingMatrix = [TrackingMatrix
                          Ellipse.toTrackingMatrix(TimeSlices{i})];
    end
    
    print('Tracking paths');
    % track with 2 slice memory, ignoring all paths with less than 3 recorded
    % positions
    Tracked = track(TrackingMatrix, 3, struct('mem', 2, 'dim', 2, ...
                                              'good', 3, 'quiet', false));
    clear TrackingMatrix
    
    % Read and merge the output of the tracking algorithm into instances of
    % EllipsePath
    print('Merging bacteria paths');
    Paths = EllipsePath.empty;
    i = 1;
    while true
        Ellipses = Ellipse.empty;
        Times = [];
        currentId = Tracked(i, 5);
        
        while i <= size(Tracked, 1) && currentId == Tracked(i, 5)
            index = Tracked(i, 3);
            t = Tracked(i, 4);
            Ellipses(end + 1) = TimeSlices{t}(index);
            Times(end + 1) = t;
            
            i = i + 1;
        end
        
        Paths(end + 1) = EllipsePath(Ellipses, Times);
        
        if i > size(Tracked, 1), break, end
    end   
end