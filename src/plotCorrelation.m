function plotCorrelation(Frames)  
    hold on;
    
    Sums = zeros(size(Frames{1}{2}));
    for i = 1:length(Frames)
        plot(Frames{i}{1}, Frames{i}{2}, 'k.', 'MarkerSize', 5);
        Sums = Sums + Frames{i}{2};
    end
    
    plot(Frames{1}{1}, Sums./length(Frames), 'b-', 'LineWidth', 3);
end