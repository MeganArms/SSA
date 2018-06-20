function [sf,t,sem,b,nume,firstframe,lastframe] = crtd(Nframes,exptime,PBflag,spots,objs)

% PB flag is true for photobleaching analysis, false otherwise
        
l = 1; r = 1;
L = cellfun(@length,spots);
reside = zeros(sum(L),1);
firstframe = reside; lastframe = reside;
for q = 1:length(spots)
    spotEvents = spots{q};
    objs_link = objs{q};
    if q == 1
        idx = r:length(spotEvents);
        r = r + length(spotEvents);
    else
        idx = r:(length(spotEvents)+r-1);
        r = r + length(spotEvents);
    end
    for k = 1:length(spotEvents)
        traj = spotEvents(k).trajectory;
        inds = find(~isnan(traj));
        t = idx(k);
        firstframe(t) = objs_link(5,traj(min(inds)));
        lastframe(t) = objs_link(5,traj(max(inds)));
        n = lastframe(t) - firstframe(t);
        if PBflag && firstframe(t) == 1 && lastframe(t) ~= Nframes && n > 0
            reside(l) = (n+1)*exptime;
            l = l + 1;
        elseif ~PBflag && firstframe(t) ~= 1 && lastframe(t) ~= Nframes && n > 0
            reside(l) = (n+1)*exptime;
            l = l + 1;
        end
    end
end
reside = reside(reside > 0);
firstframe = firstframe(reside > 0);
lastframe = lastframe(reside > 0);

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