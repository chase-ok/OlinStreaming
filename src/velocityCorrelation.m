function f = velocityCorrelation(varargin)
    function vc = calculate(Particle, TXYAdXdY)
        V = toUnit(Particle(5:6));
        vc = 0;
        
        n = size(TXYAdXdY, 1);
        for i = 1:n
            vc = vc + abs(dot(V, toUnit(TXYAdXdY(i, 5:6))));
        end
        
        if n ~= 0
            vc = vc/n;
        end
    end
    f = byRadius(@calculate, false, varargin{:});
end

function V = toUnit(V)
    m = norm(V);
    if m < 0.000001
        V = V.*0;
    else
        V = V./m;
    end
end