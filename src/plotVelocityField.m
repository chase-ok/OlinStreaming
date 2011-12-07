function h = plotVelocityField(Positions, Velocities)
% h = plotVelocityField(Positions, Velocities)
% Plots a vector field of velocities based on the output of
% calculateVelocityField.
%
% Arguments:
%     - Positions: n-by-2 matrix of positions. Column 1 is x, Column 2 is y.
%     - Velocities: n-by-2 matrix of velocities. Column 1 is Vx, Column 2 is Vy.
%
% Returns:
%     - h (optional): a handle to the resulting quiver plot object.
% 
% Usage:
%     [P V] = calculateVelocityField(imageSeries, M, 32, [1 2]);
%     plotVelocityField(P, V);

    h_ = quiver(Positions(:, 1), Positions(:, 2), ...
                Velocities(:, 1), Velocities(:, 2));
    if nargout == 1, h = h_; end
end