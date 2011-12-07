function Ellipses = findCyanoBacteria(Image, z, debug)
% Ellipses = findCyanoBacteria(Image)
% Finds the (roughly) circular cyano bacteria in a grayscale image.
%
% Arguments:
%     - Image: the cyano bacteria image (n-by-m matrix containing pixel data)
%     - z: the current z-depth to be assigned to each bacteria
%     - debug: if specified and true, the bacteria that are found are plotted on
%     top of the image for verification.
% 
% Returns:
%     - Ellipses: a vector of ellipses that were found in the image that
%     represent individual bacteria.
% 
% Usage:
%     Ellipses = findCyanoBacteria(Image, 0, true);
% 

    Filters = {@medianFilter, ...
               @otsuFilter, ...
               @medianFilter, ...
               makeSizeFilter(30, 200), ...
              };
    Filtered = logical(applyFilters(Image, Filters));
    
    Props = regionprops(Filtered, 'Centroid', 'Area', 'Orientation', ...
                        'MinorAxisLength', 'MajorAxisLength');
    Ellipses = Ellipse.fromRegionProps(Props, z);

    if nargin == 3 && debug
        plotBacteria(Image, Ellipses);
    end
end

function plotBacteria(Image, Ellipses)
    figure;
    imshow(Image);
    hold on;

    for i = 1:length(Ellipses)
        Ellipses(i).plot('b-', true);
    end
end