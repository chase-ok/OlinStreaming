function Image = applyFilters(Image, Filters)
% Image = applyFilters(Image, Filters)
% Applies a series of filters to an image and returns the final image.
%
% Arguments:
%     - Image: an image to be filtered (n-by-m matrix containing pixel data)
%     - Filters: cell vector of filtering functions. Each function should take 1
%     argument, the image, and return the filtered image.
%
% Returns:
%     - Image: the filtered image
% 
% Usage:
%     Image = applyFilters(Image, {@deconvolve, @otsuFilter});
% 

    for i = 1 : length(Filters)
        Image = Filters{i}(Image);
    end
end