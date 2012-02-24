function f = directorCorrelation(varargin)
    function dc = calculate(Particle, TXYAdXdY)
        D = toDirector(Particle(4));
        dc = 0;
        
        n = size(TXYAdXdY, 1);
        for i = 1:n
            dc = dc + abs(dot(D, toDirector(TXYAdXdY(i, 4))));
        end
        
        if n ~= 0
            dc = dc/n;
        end
    end
    f = byRadius(@calculate, false, varargin{:});
end

function D = toDirector(angle)
    D = [cos(angle), sin(angle)];
end