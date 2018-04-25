function [sf,t,sem,b,nume] = crtd(spotEvents,Nframes,exptime,objs_link)

l = 1;
reside = zeros(length(spotEvents),1);
for k = 1:length(spotEvents)
    traj = spotEvents(k).trajectory;
    inds = find(~isnan(traj));
    firstframe = objs_link(5,traj(min(inds)));
    lastframe = objs_link(5,traj(max(inds)));
    n = lastframe - firstframe;
    if firstframe ~= 1 && lastframe ~= Nframes && n > 0
        reside(l) = (n+1)*exptime;
        l = l + 1;
    end
end
reside = reside(reside > 0);

maxT = Nframes*exptime;
edges = 2*exptime-0.5*exptime:exptime:maxT+0.5*exptime-exptime;
count = histcounts(reside,edges);

T = 2*exptime:exptime:(maxT-exptime);
c = 1./(heaviside(maxT-T).*(1-(T/maxT)));

b = 1/sum(count.*c); nume = zeros(length(T),1);
for k = 1:length(T)
    nume(k) = sum(count(k+1:end).*c(k+1:end)); 
end
CRTD = nume*b;
[sf,ind] = unique(CRTD);
t = T(ind);
sem = sqrt(sf);