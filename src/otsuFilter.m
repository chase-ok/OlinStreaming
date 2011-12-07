function Image = otsuFilter(Image)
% Image = otsuFilter(Image)
% Thresholds an image based on the threshold given by the otsu algorithm.
%
% Arguments:
%     - Image: the image to be filtered
%
% Returns:
%     - Image: the filtered image
% 
% Usage:
%     Image = otsuFilter(Image);

    threshold = graythresh(Image);
    Image = im2bw(Image, threshold);
end
