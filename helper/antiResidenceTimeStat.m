% function [f,SF,T,anti_rt,anti_rtCI,count] = antiResidenceTimeStat(traj, Nframes, varargin)
function [f,x,atl,exposure,anti_rt,anti_rtCI,cens] = antiResidenceTimeStat(traj, Nframes, varargin)
% ANTIRESIDENCETIMESTAT calculates the distribution of time between landing events.
%
% INPUT:
%   - TRAJ - input trajectory data in trajectory format (output of
%   CONNECTEVENTS or EVENTLINKS)
% OUTPUT:
%   - COUNT - number of trajectories corresponding to specific anti residence
%   time.

% Get anti-residence "trajectories"
p = 1; antitrajlengths = zeros(3*length(traj),1); cens = antitrajlengths;
for i = 1:length(traj)
    if sum(isnan(traj(i).trajectory) >= 1)
        antitrajs = regionprops(isnan(traj(i).trajectory),'Area');
        numinds = length(antitrajs);
        minp = p; maxp = p + numinds - 1;
        antitrajlengths(minp:maxp) = extractfield(antitrajs,'Area');
        p = maxp + 1;
        cens(minp:maxp) = 0;
        if isnan(traj(i).trajectory(1))
            cens(minp) = 1;
        elseif isnan(traj(i).trajectory(end))
            cens(maxp) = 1;
        end
    else
        antitrajlengths(p) = Nframes - length(traj(i).trajectory);
        cens(p) = 1;
    end
end
atl = antitrajlengths(antitrajlengths ~= 0);
atl(atl < 0) = 0;
cens = cens(antitrajlengths ~= 0);

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
% 
% maxT = exposure*(Nframes);
% deltaT = exposure;
% % count = count(1:end-1)';
% T = deltaT:deltaT:maxT;
% 
% c = 1./(heaviside(maxT-[0,T]).*(1-([0,T]/maxT)));
% % c = [1,ones(1,length(T))];
% c = c(1:end-1);
% if length(c) ~= length(count)
%     n = padarray(count,[0 abs(length(c) - length(count))],0,'post');
%     SF = abs(1-cumsum(n.*c/sum(n.*c)));
%     %sem = 0.68.*SF;%n.*c/sum(n.*c);
%     sem = sqrt(SF);
% else
%     SF = abs(1-cumsum(count.*c/sum(count.*c)));
%     %sem = 0.68.*SF;%count.*c/sum(count.*c);
%     sem = sqrt(SF);
% end
% f = fit(T',SF','exp1','Weights',sem);
% anti_rt = f.a;
% ci = confint(f);
% anti_rtCI = [f.a - ci(1), ci(2) - f.a];


[f,x] = ecdf(atl*exposure,'function','survivor','censoring',cens); 
[anti_rt,anti_rtCI] = expfit(atl*exposure,'censoring',cens);

% [f,x] = ecdf(atl*exposure,'function','survivor'); 
% [anti_rt,anti_rtCI] = expfit(atl*exposure);