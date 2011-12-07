function Image = edgeReductionFilter(Image)
% Image = edgeReductionFilter(Image)
% Removes (i.e. sets to black) the edge pixels in an image. This is done in 
% order to make separating large/intersecting objects easier. Provided for use 
% with applyFilters.
%
% Arguments:
%     - Image: an image to be edge-reduced (n-by-m matrix containing pixel data)
%
% Returns:
%     - Image: the edge-reduced image
% 
% Usage:
%     Image = edgeReductionFilter(Image);
% 

    Edges = edge(Image, 'canny'); % Canny method also finds 'weak' edges, ie the
                                  % ones that we actually want to remove.
    Image(Edges) = 0;
end