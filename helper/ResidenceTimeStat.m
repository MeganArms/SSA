function [count, edges] = ResidenceTimeStat(traj,varargin)
% RESIDENCETIMESTAT calculates the distribution of residence time.
%
% INPUT:
%   - TRAJ - input trajectory data in cell format (output of
%   GETCOORDINATES)
%	- D0 - the time when the proteins were introduce to the surface, i.e. the time
%	when flow into the chamber occurred, in the format 'yyyyMMdd HHmmss'.
%	- DSTART - the time when imaging commenced in the format 'yyyyMMdd HHmmss'.
% OUTPUT:
%   - COUNT - number of trajectories corresponding to specific residence
%   time.
% t0 = datetime(d0,'InputFormat','yyyyMMdd HHmmss');
% tstart = datetime(dstart,'InputFormat','yyyyMMdd HHmmss');
% L = datenum(tstart - t0)*24*60*60;
L = 0;
reside = cellfun(@length,traj);
% Default bin edges
edges = (min(reside)-0.5):1:(max(reside)+0.5);

for pair = reshape(varargin,2,[])
    if strcmpi(pair{1},'ExposureTime') % In seconds, not milliseconds
        exposure = pair{2};
        edges = L-0.5*exposure:exposure:exposure*max(reside)+L+0.5*exposure;
    end
    if strcmpi(pair{1},'Bin')
        edges = pair{2};
    end
end

count = histcounts(reside*exposure+L,edges);