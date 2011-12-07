function f = makeLargeObjectSeparationFilter(largeThreshold, n)
% f = maeLargeObjectSeparationFilter(Window)
% Creates a function that attempts to break apart objects larger than
% largeThreshold. If a large object has a narrow enough section, it will be
% separated into 2 smaller objects.
%
% Arguments:
%     - largeThreshold: the size, in pixels, above which this function will try
%     to break apart objects.
%     - n: the number of times to apply this filter (larger numbers may be more
%     successful, but are also more likely to distort images).
%
% Returns:
%     - f: the filter function that takes 1 argument, an image, and returns the
%     filtered image.
% 
% Usage:
%     f = makeLargeObjectSeparationFilter(30, 2);
%     Image = f(Image);

    function Image = filter(Image)
        sizeFilter = makeSizeFilter(largeThreshold, largeThreshold*100);
        LargeObjects = logical(sizeFilter(Image));

        Image(LargeObjects) = 0;
        for i = 1 : n
            LargeObjects = removePerimeter(LargeObjects);
        end
        Image(LargeObjects) = 1;
    end
    
    f = @filter;
end

function Image = removePerimeter(Image)    
    shape = strel('disk', 1, 0);
    Image = imerode(Image, shape);
    %Image = bwmorph(Image,'skel',Inf); %TODO: skel won't break apart objects
    Image = imdilate(Image, shape);
end