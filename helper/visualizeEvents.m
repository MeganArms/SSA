spotEvents = connectEvents_rp(R,C,F,2);
spots = struct2table(spotEvents);
events = spots.trajectory;

% Find any indices that are visited by multiple molecules
coords = C(cellfun(@length,C)>2);
intCoords = cellfun(@(x) round(x'),coords,'UniformOutput',0);
intInds = cellfun(@(x) sub2ind([512 512],x(:,1),x(:,2)),intCoords,'UniformOutput',0);
% If the trajectory is on multiple pixels, the pixels are counted individually
uniqueInds_inTraj = cell(length(intInds),1);
indLengths = cell(length(intInds),1);
for k = 1:length(intInds)
    uniqueInds_inTraj{k} = unique(intInds{k});
    indLengths{k} = repmat(length(intInds{k}),[length(uniqueInds_inTraj{k}) 1]);
end
uniqueInds_inTraj = cell2mat(uniqueInds_inTraj); 
indLengths = cell2mat(indLengths);
% Only take the pixels that are landed on more than once
[sortedUniqueInds_inTraj, inds] = sort(uniqueInds_inTraj);
sortedIndLengths = indLengths(inds);
repeatInds_allTraj = sortedUniqueInds_inTraj(diff(sortedUniqueInds_inTraj) == 0); 
indLengths_allTraj = sortedIndLengths(diff(sortedUniqueInds_inTraj) == 0);

% Get the frequency image of locations
edges = 0.5:1:(512)^2+0.5;
[counts,edges] = histcounts(repeatInds_allTraj,edges);
counts(counts>0) = counts(counts>0)+1;
countsIm = reshape(counts,[512 512]);
% This only works because the length data is sorted in the order of the
% frequency data and the bin size is - actually the bin size doesn't matter
% lifetimefreq = zeros(length(counts)-1,1); l = 1;
% for k = 1:length(counts)-1
%     lifetimestotake = counts(k);
%     if lifetimestotake ~= 0
%         lifetimefreq(k) = mean(indLengths_allTraj(l:l+lifetimestotake-1));
%         l = l+lifetimestotake;
%     end
% end

% Visualize
zerospots = find(countsIm == 0);
clearBG = countsIm;
clearBG(zerospots) = NaN;
figure,mesh(clearBG);
xlabel('X','FontSize',14), ylabel('Y','FontSize',14)
zlabel('f(Events)','FontSize',14)
title('Event Frequency','FontSize',16);
h = gca; h.YLim = [0 512]; h.XLim = [0 512]; clear h