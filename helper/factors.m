function [icounts,dtime,firstframe,lastframe,brightness,avgb,rangeb,maxb,rcdfb] = factors(Nframes,PBflag,spots,objs)

% PB flag is true for photobleaching analysis, false otherwise

r = 1;
L = cellfun(@length,spots);
firstframe = zeros(sum(L),1);
s = 1; 
p = 1; 
lastframe = firstframe;
avgb = firstframe;
rangeb = firstframe;
dtime = firstframe;
maxb = 0;
icounts = zeros(Nframes*sum(L),4);
brightness = cell(sum(L),1);
bmat = NaN(sum(L),Nframes);
for u = 1:length(spots)
    spotEvents = spots{u};
    objs_link = objs{u};
    if u == 1
        idx = r:length(spotEvents);
        r = r + length(spotEvents);
    else
        idx = r:(length(spotEvents)+r-1);
        r = r + length(spotEvents);
    end
    for k = 1:length(spotEvents)
        traj = spotEvents(k).trajectory;
        t = idx(k);
        inds = find(~isnan(traj));
        frames = objs_link(5,traj(inds));
        int_counts = objs_link(3,traj(inds));
        bmat(t,frames) = int_counts;
        brightness{t} = cumsum(bmat(t,:),'reverse','omitnan');
        maxb = max([maxb,max(int_counts)]);
        avgb(t) = mean(int_counts);
        rangeb(t) = range(int_counts);
        firstframe(t) = frames(1);
        lastframe(t) = frames(end);
        n = lastframe(t) - firstframe(t) + 1;

        brights = ~isnan(traj);
        darks = ~brights;
        bCC = bwconncomp(brights);
        dCC = bwconncomp(darks);

        icounts_hold = struct2cell(regionprops(bCC,traj,'PixelValues'));

        if PBflag && firstframe(t) == 1 && lastframe(t) ~= Nframes && n > 0
            % If the bright events do not overlap the start or end of the movie
            % The dark events cannot be at the start or end of a trajectory.
            lengths = cellfun(@length,dCC.PixelIdxList);
            lengths = lengths(lengths>1);
            dtime(p:p+length(lengths)-1) = lengths;
            p = p + length(lengths);

            for q = 1:length(icounts_hold)
                % Total integrated counts of the bright event. Error will be sqrt.
                icounts(s,1) = sum(objs_link(3,icounts_hold{q}));
                % Length of bright event
                icounts(s,2) = length(icounts_hold{q});
                s = s + 1;
            end
            % Number of bright events in this traj in column 3
            icounts(s-q:s-1,3) = q*ones(q,1);
            % Length of this entire connected traj in column 4
            icounts(s-q:s-1,4) = n*ones(q,1);
        elseif ~PBflag && firstframe(t) ~= 1 && lastframe(t) ~= Nframes && n > 0
            % If the bright events do not overlap the start or end of the movie
            % The dark events cannot be at the start or end of a trajectory.
            lengths = cellfun(@length,dCC.PixelIdxList);
            lengths = lengths(lengths>1);
            dtime(p:p+length(lengths)-1) = lengths;
            p = p + length(lengths);

            for q = 1:length(icounts_hold)
                % Total integrated counts of the bright event. Error will be sqrt.
                icounts(s,1) = sum(objs_link(3,icounts_hold{q}));
                % Length of bright event
                icounts(s,2) = length(icounts_hold{q});
                s = s + 1;
            end
            % Number of bright events in this traj in column 3
            icounts(s-q:s-1,3) = q*ones(q,1);
            % Length of this entire connected traj in column 4
            icounts(s-q:s-1,4) = n*ones(q,1);
        end
    end
end
avgAllB = mean(bmat,1,'omitnan');
rcdfb = cumsum(avgAllB,'reverse','omitnan');
icounts = icounts(1:s-1,:);
dtime = dtime(1:p-1);