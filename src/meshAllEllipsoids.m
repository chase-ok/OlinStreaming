function meshAllEllipsoids(Ellipsoids)
% meshAllEllipsoids(Ellipsoids)
% Takes a vector of Ellipsoids and creates a 3D plot.
%
% Arguments:
%     - Ellipsoids: vector of Ellipsoid
% 
% Usage:
%     meshAllEllipsoids(Ellipsoids)

    figure; hold on;
    for i = 1 : length(Ellipsoids)
        meshEllipses(Ellipsoids(i));
    end
    axis([0 512 0 512 0 max([Ellipsoids.z])]);
    xlabel('x'); ylabel('y'); zlabel('z');
end