function Ellipses = findEcoliBacteria(Image, z, debug)
% Ellipses = findEcoliBacteria(Image)
% Finds the (roughly) ellipse E. coli bacteria in a grayscale image.
%
% Arguments:
%     - Image: the E. coli bacteria image (n-by-m matrix containing pixel data)
%     - z: the current z-depth to be assigned to each bacteria
%     - debug: if specified and true, the bacteria that are found are plotted on
%     top of the image for verification.
% 
% Returns:
%     - Ellipses: a vector of ellipses that were found in the image that
%     represent individual bacteria.
% 
% Usage:
%     Ellipses = findEcoliBacteria(Image, 0, true);
% 
    %removeBackground, 
    Filters = {@removeBackground, ...
               makeThresholdFilter(0.18), ... % very sensitive, need to fix
               makeSizeFilter(4, 100), ...
              };
    
%     Filters = {@deconvolve, ...
%                @edgeReductionFilter, ...
%                makeThresholdFilter(0.8), ... 
%                makeSizeFilter(1.5, 200), ...
%               };
    Filtered = logical(applyFilters(Image, Filters));
    
    Props = regionprops(Filtered, 'Centroid', 'Area', 'Orientation', ...
                        'MinorAxisLength', 'MajorAxisLength');
    Ellipses = Ellipse.fromRegionProps(Props, z);
    %Ellipses = Ellipse.splitLargeEllipses(Ellipses);
    
    if nargin == 3 && debug
        plotBacteria(Image, Ellipses);
    end
end

function plotBacteria(Image, Ellipses)
    figure;
    imshow(imadjust(Image));
    hold on;

    for i = 1:length(Ellipses)
        Ellipses(i).plot('b-', true);
    end
end