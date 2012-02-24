function f = particleDistance(varargin)
    f = byRadius(@(p, Particles) size(Particles, 1), true, varargin{:});

%     if nargin == 0
%         Radiuses = 1:3:50;
%     elseif nargin == 1
%         Radiuses = dr;
%     elseif nargin == 2
%         Radiuses = 0.5:dr:maxRadius;
%     end
% 
%     RadiusesSq = Radiuses.^2;
%     Areas = pi.*[RadiusesSq(1), (RadiusesSq(2:end) - RadiusesSq(1:end-1))];
%     
%     function Result = calculate(imageSeries, TXYAdXdY)
%         XY = TXYAdXdY(:, 2:3);
%         n = size(XY, 1);
%         
%         NumParticles = zeros(size(Radiuses));
%         
%         for i = 1:n
%             DistancesSq = sum((XY - repmat(XY(i, :), n, 1)).^2, 2);
%             
%             % Use accum to keep track of the number of particles within the 
%             % previous radius so we can just subtract it off instead of doing
%             % two distance comparisons per particle.
%             accum = 1; % 1 for the particle we're measuring from
%             for j = 1:(length(RadiusesSq) - 1)
%                 rSq = RadiusesSq(j);
%                 numInside = sum(DistancesSq <= rSq) - accum;
%                 NumParticles(j) = NumParticles(j) + numInside;
%                 accum = accum + numInside;
%             end
%         end
%         
%         % Take care of the remaining particles
%         %Radiuses(end + 1) = Radiuses(end) + 1;
%         %NumParticles(end + 1) = n - accum;
%         
%         % Dividing by the total number of particles (size(XY, 1)) gives the mean
%         % number of particles within a radius, dividing again gives the 
%         % probability that a given particle is that far away.
%         Probability = NumParticles./(accum^2);
%         Density = Probability./Areas;
%         
%         Result = {Radiuses, Density};
%     end
%     
%     f = @calculate;
end