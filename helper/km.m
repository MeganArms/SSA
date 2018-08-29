function [f,x,reside,cens,flo,fup] = km(Nframes,exptime,PBflag,spots)

% PB flag is true for photobleaching analysis, false otherwise

l = 1; r = 1;
L = cellfun(@length,spots);
reside = zeros(sum(L),1);
cens = reside;
firstframe = reside; lastframe = reside;
for q = 1:length(spots)
    spotEvents = spots{q};
    if q == 1
        idx = r:length(spotEvents);
        r = r + length(spotEvents);
    else
        idx = r:(length(spotEvents)+r-1);
        r = r + length(spotEvents);
    end
    for k = 1:length(spotEvents)
        frames = spotEvents(k).frames;
        t = idx(k);
        firstframe(t) = min(frames);
        lastframe(t) = max(frames);
        n = lastframe(t) - firstframe(t);
        if PBflag && firstframe(t) == 1 && n > 0
            reside(l) = (n+1)*exptime;
            if lastframe(t) == Nframes
                cens(l) = 1;
            else
                cens(l) = 0;
            end
            l = l + 1;
        elseif ~PBflag && firstframe(t) ~= 1 && n > 0
            reside(l) = (n+1)*exptime;
            if lastframe(t) == Nframes
                cens(l) = 1;
            else
                cens(l) = 0;
            end
            l = l + 1;
%         elseif ~PBflag && n > 0
%             reside(l) = (n+1)*exptime;
%             if lastframe(t) == Nframes || firstframe(t) == 1
%                 cens(l) = 1;
%             else
%                 cens(l) = 0;
%             end
%             l = l + 1;
        end
    end
end
reside = reside(reside > 0);
cens = cens(reside > 0);

[f,x,flo,fup] = ecdf(reside,'function','survivor','censoring',cens);