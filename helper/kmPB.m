function [f,x,reside,cens] = kmPB(spotEvents1,Nframes,exptime,objs_link1)

l = 1;
reside = zeros(length(spotEvents1),1);
cens = reside;
firstframe = reside; lastframe = reside;
for k = 1:length(spotEvents1)
    traj = spotEvents1(k).trajectory;
    inds = find(~isnan(traj));
    firstframe(k) = objs_link1(5,traj(min(inds)));
    lastframe(k) = objs_link1(5,traj(max(inds)));
    n = lastframe(k) - firstframe(k);
    if firstframe(k) == 1 && n > 0
        reside(l) = (n+1)*exptime;
        if lastframe(k) == Nframes
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