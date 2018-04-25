function [f,x,reside,cens] = km2(spotEvents1,spotEvents2,Nframes,exptime,objs_link1,objs_link2)

l = 1;
reside = zeros(length(spotEvents1)+length(spotEvents2),1);
cens = reside;
firstframe = reside; lastframe = reside;
for k = 1:length(spotEvents1)
    traj = spotEvents1(k).trajectory;
    inds = find(~isnan(traj));
    firstframe(k) = objs_link1(5,traj(min(inds)));
    lastframe(k) = objs_link1(5,traj(max(inds)));
    n = lastframe(k) - firstframe(k);
    if firstframe(k) ~= 1 && n > 0
        reside(l) = (n+1)*exptime;
        if lastframe(k) == Nframes
            cens(l) = 1;
        else
            cens(l) = 0;
        end
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
    if firstframe(t) ~= 1 && n > 0
        reside(l) = (n+1)*exptime;
        if lastframe(t) == Nframes
            cens(l) = 1;
        else
            cens(l) = 0;
        end
        l = l + 1;
    end
end
reside = reside(reside > 0);
cens = cens(reside > 0);

[f,x] = ecdf(reside,'function','survivor','censoring',cens);