function Image = invertFilter(Image)
% Image = invertFilter(Image)
% Inverts an image (the brightest pixel becomes the darkest and vice versa). 
% Provided for use with applyFilters.
%
% Arguments:
%     - Image: an image to be inverted (n-by-m matrix containing pixel data)
%
% Returns:
%     - Image: the inverted image
% 
% Usage:
%     Image = invertFilter(Image);

    Image = max(max(Image)) - Image;
end