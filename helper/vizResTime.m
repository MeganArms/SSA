% Get indices of all trajectories with length > 2 frames
longEvents = C(cellfun(@length,C) > 2);

% Get the lifetime for each pixel with an even on it
y = cellfun(@(x) round(x'),longEvents,'UniformOutput',0);
z = cellfun(@(x) sub2ind([512 512],x(:,1),x(:,2)),y,'UniformOutput',0);
inds = cellfun(@unique,z,'UniformOutput',0);
LTs = cellfun(@(x) repmat(length(x)*2,length(unique(x)),1),z,'UniformOutput',0);
lifetimes = cell2mat([inds, LTs]);

uniqueInds = unique(lifetimes(:,1));
indVect = zeros(512^2,1);
for k = 1:length(uniqueInds)
    indVect(uniqueInds(k)) = mean(lifetimes(lifetimes(:,1)==uniqueInds(k),2));
end
ltIm = reshape(indVect,[512 512]);

