function [icounts,dtime,firstframe,lastframe,brightness,avgb,rangeb,maxb,rcdfb] = factorsPB(spotEvents1,Nframes,objs_link1)


% Including only the bright and dark times that do not overlap the start or
% finish of the video, which is consistent with the KS correction.
s = 1; 
p = 1; 
% l = 1;
firstframe = zeros(length(spotEvents1),1);
lastframe = firstframe;
avgb = firstframe;
rangeb = firstframe;
maxb = 0;
icounts = zeros(Nframes*(length(spotEvents1)),4);
dtime = zeros(Nframes*(length(spotEvents1)),1);
brightness = cell(length(spotEvents1),1);
bmat = NaN(length(spotEvents1),Nframes);
for k = 1:length(spotEvents1)
    traj = spotEvents1(k).trajectory;
    if eventfreq1(k) > 1
        continue
    else
    inds = find(~isnan(traj));
    frames = objs_link1(5,traj(inds));
    int_counts = objs_link1(3,traj(inds));
    bmat(k,frames) = int_counts;
%     brightness{k} = cumsum(bmat(k,:),'reverse','omitnan');
    brightness{k} = bmat(k,:);
    maxb = max([maxb,max(int_counts)]);
    avgb(k) = mean(int_counts);
    rangeb(k) = range(int_counts);
    firstframe(k) = frames(1);
    lastframe(k) = frames(end);
    n = lastframe(k) - firstframe(k) + 1;
    
    brights = ~isnan(traj);
    darks = ~brights;
    bCC = bwconncomp(brights);
    dCC = bwconncomp(darks);
     
    icounts_hold = struct2cell(regionprops(bCC,traj,'PixelValues'));
    
    if firstframe(k) == 1 && lastframe(k) ~= Nframes && n > 0
        % If the bright events do not overlap the end of the movie
%         btime(l:l+length(bCC.PixelIdxList)-1) = cellfun(@length,bCC.PixelIdxList);
        lengths = cellfun(@length,dCC.PixelIdxList);
        dtime(p:p+length(lengths)-3) = lengths(2:end-1);
%         l = l + length(bCC.PixelIdxList);
        p = p + length(lengths) - 2;
        
        for q = 1:length(icounts_hold)
            % Total integrated counts of the bright event. Error will be sqrt.
            icounts(s,1) = sum(objs_link1(3,icounts_hold{q}));
            % Length of bright event
            icounts(s,2) = length(icounts_hold{q});
            s = s + 1;
        end
        % Number of bright events in this traj in column 3
        icounts(s-q:s-1,3) = q*ones(q,1);
        % Length of this entire connected traj in column 4
        icounts(s-q:s-1,4) = n*ones(q,1);
%     elseif firstframe(k) ~= 1 && lastframe(k) == Nframes && n > 1
%         % If the last bright event appears on the last frame
% %         blengths = cellfun(@length,bCC.PixelIdxList);
% %         btime(l:l+length(blengths)-2) = blengths(1:end-1);
%         dlengths = cellfun(@length,dCC.PixelIdxList);
%         dtime(p:p+length(dlengths)-2) = dlengths(2:end);
% %         l = l + length(blengths) - 1;
%         p = p + length(dlengths) - 1;
%         
%         if length(icounts_hold) > 1
%             for q = 1:length(icounts_hold)-1
%                 icounts(s,1) = sum(objs_link1(3,icounts_hold{q}));
%                 icounts(s,2) = length(icounts_hold{q});
%                 s = s + 1;
%             end
%             icounts(s-q:s-1,2) = q*ones(q,1);
%             icounts(s-q:s-1,3) = n*ones(q,1);
%         end
%     elseif firstframe(k) == 1 && lastframe(k) ~= Nframes && n > 1
%         % If the first bright event is on the first frame
% %         blengths = cellfun(@length,bCC.PixelIdxList);
% %         btime(l:l+length(blengths)-2) = blengths(2:end);
%         dlengths = cellfun(@length,dCC.PixelIdxList);
%         dtime(p:p+length(dlengths)-2) = dlengths(1:end-1);
% %         l = l + length(blengths) - 1;
%         p = p + length(dlengths) - 1;
%         
%         if length(icounts_hold) > 1
%             for q = 2:length(icounts_hold)
%                 icounts(s,1) = sum(objs_link1(3,icounts_hold{q}));
%                 icounts(s,2) = length(icounts_hold{q});
%                 s = s + 1;
%             end
%             icounts(s-q:s-1,2) = (q-1)*ones(q-1,1);
%             icounts(s-q:s-1,3) = n*ones(q,1);
%         end
%     elseif firstframe(k) == 1 && lastframe(k) == Nframes && n > 1
%         % If the bright events overlap first and last frames
% %         blengths = cellfun(@length,bCC.PixelIdxList);
% %         btime(l:l+length(blengths)-2) = blengths(2:end);
%         dlengths = cellfun(@length,dCC.PixelIdxList);
%         dtime(p:p+length(dlengths)-1) = dlengths(1:end);
% %         l = l + length(blengths) - 1;
%         p = p + length(dlengths) - 1;
%         
%         if length(icounts_hold) > 2
%             for q = 2:length(icounts_hold)-1
%                 icounts(s,1) = sum(objs_link1(3,icounts_hold{q}));
%                 icounts(s,2) = length(icounts_hold{q});
%                 s = s + 1;
%             end
%             icounts(s-q:s-1,2) = (q-1)*ones(q-1,1);
%             icounts(s-q:s-1,3) = n*ones(q,1);
%         end
%     else
%         disp('something is wrong')
    
    end
    end
end

avgAllB = mean(bmat,1,'omitnan');
rcdfb = cumsum(avgAllB,'reverse','omitnan');
icounts = icounts(1:s-1,:);
dtime = dtime(1:p-1);
