function Image = deconvolve(Image)
% Image = deconvolve(Image)
% Deconvolves the image using Lucy-Richardson method and settings appropriate to
% Olin's confocal. Provided for use with applyFilters.
%
% Arguments:
%     - Image: an image to be deconvolved (n-by-m matrix containing pixel data)
%
% Returns:
%     - Image: the deconvolved image
% 
% Usage:
%     Image = deconvolve(Image);
% 
    Image = deconvlucy(Image, fspecial('gaussian', 4, 4));
end