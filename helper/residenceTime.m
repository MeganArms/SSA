function [f,x,SF,count,sem,CC,b,xx,nume] = residenceTime(traj,frames,Nframes,exposure)
% RESIDENCETIMESTAT calculates the distribution of residence time.
%
% INPUT:
%   - TRAJ - input trajectory data in cell format (output of
%   COLLATE)
% OUTPUT:
%   - COUNT - number of trajectories corresponding to specific residence
%   time.
toexclude = cellfun(@(x)x(1)==1||x(end)==Nframes,frames);
newtraj = traj(~toexclude);
lengths = cellfun(@length,newtraj);
reside = lengths(lengths>1);
% Get counts
maxT = exposure*(Nframes);
edges = 2*exposure-0.5*exposure:exposure:maxT+0.5*exposure-exposure;
count = histcounts(reside*exposure,edges);

T = 2*exposure:exposure:(maxT-exposure);
c = 1./(heaviside(maxT-T).*(1-(T/maxT)));
SF = 1-cumsum(count.*c/sum(count.*c));
sem = sqrt(count);

% Test another method
b = 1/sum(count.*c); nume = zeros(length(T),1);
for k = 1:length(T)
    nume(k) = sum(count(k+1:end).*c(k+1:end)); 
end
CRTD = nume*b;
[CC,ind] = unique(CRTD);
xx = T(ind);
sem = sem(ind);

% Compare to KM
toexclude = cellfun(@(x)x(1)==1,frames);
f = frames(~toexclude);
newtraj = traj(~toexclude);
lengths = cellfun(@length,newtraj);
reside = lengths(lengths>1);
cens = cellfun(@(x)x(end)==Nframes,f);
[f,x] = ecdf(reside*exposure,'function','survivor','censoring',cens);