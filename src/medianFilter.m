function Image = medianFilter(Image)
% Image = medianFilter(Image)
% Applies a median filter to the image.
%
% Arguments:
%     - Image: the image to be filtered
%
% Returns:
%     - Image: the filtered image
% 
% Usage:
%     Image = medianFilter(Image);

    Image = medfilt2(Image); % oh hey look, didn't need to implement it
                             % ourselves after all...
end