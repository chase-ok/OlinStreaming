function f = makeThresholdFilter(threshold)
% f = makeThresholdFilter(threshold)
% Creates a function that thresholds a grayscale image into a binary image.
%
% Arguments:
%     - threshold: a pixel value (0.0 to 1.0) above which pixels are rendered 
%     as 1, below which they are rendered as 0.
%
% Returns:
%     - f: the filter function that takes 1 argument, an image, and returns the
%     filtered image.
% 
% Usage:
%     f = makeThresholdFilter(0.7);
%     Image = f(Image);

    f = @(Image) im2bw(Image, threshold);
end