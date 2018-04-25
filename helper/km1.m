function [f,x,reside,cens] = km1(spotEvents,Nframes,exptime)

l = 1;
reside = zeros(length(spotEvents),1);
cens = reside;
% firstframe = reside; lastframe = reside;
for k = 1:length(spotEvents)
%     traj = spotEvents(k).trajectory;
    frms = spotEvents(k).frames;
%     
%     inds = find(~isnan(traj));
%     firstframe(k) = objs_link(5,traj(min(inds)));
%     lastframe(k) = objs_link(5,traj(max(inds)));
    n = frms(end) - frms(1);
    if frms(1) ~= 1 && n > 0
        reside(l) = (n+1)*exptime;
        if frms(end) == Nframes
            cens(l) = 1;
        else
            cens(l) = 0;
        end
        l = l + 1;
    end
end
cens = cens(reside > 0);
reside = reside(reside > 0);


[f,x] = ecdf(reside,'function','survivor','censoring',cens);