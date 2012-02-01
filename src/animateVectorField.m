function animateVectorField(imageSeries, Results, pauseDuration)
    if nargin < 3, pauseDuration = 0.05; end

    for frame = 1:length(Results)
        P = Results{frame}{1};
        V = Results{frame}{2};
        
        quiver(P(:, 1), P(:, 2), V(:, 1), V(:, 2));
        axis([0 imageSeries.ImageSize(1), 0 imageSeries.ImageSize(2)]);
        drawnow;
        pause(pauseDuration);
    end
    
end