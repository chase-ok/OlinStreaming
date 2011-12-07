function Images = loadImageSeries(first, lastImage, last)
% Image = loadImageSeries(first, lastImage, last)
% Utility function to make loading numbered image lists easier. Loads series
% where the file names are formated like [first]NNN[last] where N are digits
% (2, 3, or 4 digits)
%
% Arguments:
%     - first: the first part of the file name of every image in the series as a
%     string.
%     - lastImage: integer, the number of the last image in the series.
%     - last: the last part of the file name of every image in the series as a
%     string.
%
% Returns:
%     - Images: a cell array of images (n-by-n matrices of pixel value).
% 
% Usage:
%     Images = loadImageSeries('Series_t', 80, '_ch00.tif');

    if lastImage < 100
        formatString = '%02u';
    elseif lastImage < 1000
        formatString = '%03u';
    else
        formatString = '%04u';
    end
    
    for i = 0 : lastImage
        i
        Images{i + 1} = imread([first, num2str(i, formatString), last]);
    end
end