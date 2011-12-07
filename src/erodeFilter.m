function Image = erodeFilter(Image)
% Image = erodeFilter(Image)
% Morphological erosion with a diamond shaped element. Provided for use with 
% applyFilters.
%
% Arguments:
%     - Image: an image to be eroded (n-by-m matrix containing pixel data)
%
% Returns:
%     - Image: the eroded image
% 
% Usage:
%     Image = erodeFilter(Image);
% 

    Image = ordfilt2(Image, 2, [0 1 0
                                1 1 1
                                0 1 0]);
end