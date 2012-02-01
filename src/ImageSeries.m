classdef ImageSeries < handle
    
    properties
        ImageFiles % Cell array of file paths
        Images % Cell array of images
        Pixel % [x, y] pixel size in microns
        ImageSize % in pixels
        dt % in seconds
        startTime % datevec
        seriesNum
        channelNum
        folder
        numFrames
    end
    
    properties (Constant)
        TRACKING_PY_PATH = '/home/ckernan/dev/OlinStreaming/src/python/tracking.py';
    end
    
    methods
        function obj = ImageSeries(folder, seriesNum, channelNum)
            if nargin == 0, return, end
            
            obj.folder = folder;
            obj.seriesNum = seriesNum;
            obj.channelNum = channelNum;
        end
        
        function load(obj)
        % Parses the properties file and loads the images into the file.
        
            obj.parseTrackingPyOutput();
            obj.loadImages();
            obj.numFrames = length(obj.Images);
            obj.ImageSize = size(obj.Images{1}).*obj.Pixel;
        end
        
        function releaseImages(obj)
        % Free up memory by clearing the Images variables.
        
            obj.Images = {};
            obj.ImageFiles = {};
        end
            
    end
        
    methods (Hidden)
        
        function parseTrackingPyOutput(obj)
        % Use the tracking.py script to parse basic information about a series,
        % and extract a list of  image file names.
        
            [status, dataStr] = system(['python "' obj.TRACKING_PY_PATH '" "', ...
                                        obj.folder '" ', ...
                                        int2str(obj.seriesNum) ' ', ...
                                        int2str(obj.channelNum)
                                       ]);
            if status ~= 0
                error(['tracking.py failed: ', dataStr]);
            end

            Lines = regexp(dataStr, '\n', 'split');
            
            obj.startTime = datevec(Lines(1));
            obj.dt = str2double(Lines(2));
            
            obj.Pixel(1) = str2double(Lines(3));
            obj.Pixel(2) = str2double(Lines(4));

            % The result of the lines (minus the last blank one) are image files
            %(in order)
            obj.ImageFiles = Lines(5:end-1); % pull off the last blank str
        end
        
        function loadImages(obj)
            obj.Images = cell(1, length(obj.ImageFiles));
            for j = 1 : length(obj.ImageFiles)
                obj.Images{j} = imread(obj.ImageFiles{j});
            end
        end
            
    end
    
end

