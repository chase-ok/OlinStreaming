function f = makeContrastFilter(Window)
% f = makeContrastFilter(Window)
% Creates a function that filters an image to only include pixels within the
% given window. Pixels with a value beneath the lower limit of the window become
% 0, pixels above the upper limit become 1, and anything between is scaled
% between the window edges.
%
% Arguments:
%     - Window: a vector with two elements: (1) lower pixel limit, (2) upper
%     pixel limit.
%
% Returns:
%     - f: the filter function that takes 1 argument, an image, and returns the
%     filtered image.
% 
% Usage:
%     f = makeContrastFilter([0.3 0.7]);
%     Image = f(Image);

    function Image = contrast(Image)
        Image = imadjust(Image, Window, []);
    end
    f = @contrast;
end