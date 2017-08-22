function [count, edges, atl] = antiResidenceTimeStat(traj,varargin)
% ANTIRESIDENCETIMESTAT calculates the distribution of time between landing events.
%
% INPUT:
%   - TRAJ - input trajectory data in trajectory format (output of
%   CONNECTEVENTS)
% OUTPUT:
%   - COUNT - number of trajectories corresponding to specific anti residence
%   time.

% Get anti-residence "trajectories"
p = 1; i = 1; antitrajlengths = zeros(3*length(traj),1);
while i < length(traj)
    if sum(isnan(traj(i).trajectory) >= 1)
        antitrajs = regionprops(isnan(traj(i).trajectory),'Area');
        numinds = length(antitrajs);
        minp = p; maxp = p + numinds - 1;
        antitrajlengths(minp:maxp) = extractfield(antitrajs,'Area');
        p = maxp + 1;
    end
    i = i + 1;
end
atl = antitrajlengths(antitrajlengths ~= 0);

% Default bin edges
edges = (min(atl)-0.5):1:(max(atl)+0.5);

for pair = reshape(varargin,2,[])
    if strcmpi(pair{1},'ExposureTime') % In seconds, not milliseconds
        exposure = pair{2};
        edges = 2*exposure-0.5*exposure:exposure:exposure*max(atl)+0.5*exposure;
    end
    if strcmpi(pair{1},'Bin')
        edges = pair{2};
    end
end

count = histcounts(atl*exposure,edges);