classdef Claim < handle
% Class used by ellipsoid detection to represent many ellipsoid's claim on an
% ellipse.

    properties
        object % the object be claimed
        Stakes % structure with 3 cell arrays: claimer, score, rank
    end
    
    properties (Dependent)
        hasStakes % true iff at least 1 object has a stake in this claim
    end
    
    methods
        function obj = Claim(object)
            if nargin == 0, return, end
            
            obj.object = object;
            obj.Stakes = struct('claimer', {}, 'score', {}, 'rank', {});
        end
        
        function has = get.hasStakes(obj)
            has = ~isempty(obj.Stakes);
        end
        
        function obj = stake(obj, claimer, score, rank)
        % Adds a stake to this claim. Claimer is the object holding the stake,
        % score represents how good of a fit this claim is (higher is better),
        % and rank is how this claim ranks against the claimer's other ranks (1
        % means this is the claimer's top claim, etc.).
            obj.Stakes(end + 1) = struct('claimer', claimer, 'score', score, ...
                                         'rank', rank);
        end
        
        function revoke(obj, claimer)
        % Remove a claimer from the stakes.
            obj.Stakes([obj.Stakes.claimer] == claimer) = [];
        end
        
        function claimer = getWinner(obj)
        % Picks the claimer best suited for this claim. Sorts first by rank,
        % then by score.
            if (~obj.hasStakes), error('No stakes'), end
                
            Points = [obj.Stakes.rank; obj.Stakes.score]';
            [~, Indices] = sortrows(Points, [1, -2]); % sort rank then score
            claimer = obj.Stakes(Indices(1)).claimer;
        end
    end
    
end

