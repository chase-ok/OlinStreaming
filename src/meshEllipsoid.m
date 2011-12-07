function meshHandle = meshEllipsoid(e)
% meshHandle = meshEllipsoid(ellipsoid)
% Creates a 3D mesh of an ellipsoid.
%
% Arguments:
%     - e: an ellipsoid
%
% Returns:
%     - meshHandle: a handle to a mesh object.
% 
% Usage:
%     h = meshEllipsoid(ellipsoid)

    [X Y Z] = ellipsoid(e.x, e.y, e.z, e.a, e.b, e.b, 10);
    meshHandle = mesh(X, Y, Z);
    rotate(meshHandle, [1 0 0], e.phi*180/pi);
    rotate(meshHandle, [0 1 0], e.theta*180/pi - 90);
end