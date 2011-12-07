function f = makeSizeFilter(minSize, maxSize)
% f = makeSizeFilter(Window)
% Creates a function that filters an image to only objects that have a pixel
% area between minSize and maxSize.
%
% Arguments:
%     - minSize: minimum object size in pixels
%     - maxSize: maximum object size in pixels
%
% Returns:
%     - f: the filter function that takes 1 argument, an image, and returns the
%     filtered image.
% 
% Usage:
%     f = makeSizeFilter(5, 30);
%     Image = f(Image);

    function Image = filter(Image)
        PixelObjects = bwconncomp(Image);
        Image = zeros(size(Image));

        for j = 1 : PixelObjects.NumObjects
            objSize = length(PixelObjects.PixelIdxList{j});
            if objSize > minSize && objSize < maxSize
                Image(PixelObjects.PixelIdxList{j}) = 1;
            end
        end
    end

    f = @filter;
end