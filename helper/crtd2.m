function [sf,t,sem,b,nume,firstframe,lastframe] = crtd2(spotEvents1,spotEvents2,Nframes,exptime,objs_link1,objs_link2)

l = 1;
reside = zeros(length(spotEvents1)+length(spotEvents2),1);
firstframe = reside; lastframe = reside;
for k = 1:length(spotEvents1)
    traj = spotEvents1(k).trajectory;
    inds = find(~isnan(traj));
    firstframe(k) = objs_link1(5,traj(min(inds)));
    lastframe(k) = objs_link1(5,traj(max(inds)));
    n = lastframe(k) - firstframe(k);
    if firstframe(k) ~= 1 && lastframe(k) ~= Nframes && n > 0
        reside(l) = (n+1)*exptime;
        l = l + 1;
    end
end
idx = 1+length(spotEvents1):(length(spotEvents1)+length(spotEvents2));
for k = 1:length(spotEvents2)
    traj = spotEvents2(k).trajectory;
    t = idx(k);
    inds = find(~isnan(traj));
    firstframe(t) = objs_link2(5,traj(min(inds)));
    lastframe(t) = objs_link2(5,traj(max(inds)));
    n = lastframe(t) - firstframe(t);
    if firstframe(t) ~= 1 && lastframe(t) ~= Nframes && n > 0
        reside(l) = (n+1)*exptime;
        l = l + 1;
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