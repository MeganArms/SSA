
% Get indices of all trajectories with length > 2 frames
coords = cell(length(events),1);
for k = 1:length(events)
    coords{k}(1:2,:) = objs_link(1:2,events{k}(~isnan(events{k})));
end

longEvents = coords(cellfun(@length,coords) > 2);
longEvents = cellfun(@(x) x',longEvents,'UniformOutput',0);
% Get residence time on each spot

% Convert to indices
allLocs = round(cell2mat(longEvents));
allInds = sub2ind([512 512], allLocs(:,1),allLocs(:,2));

% Get the frequency image of locations
edges = 0.5:1:(512)^2+0.5;
[counts,edges] = histcounts(allInds,edges);
countsIm = reshape(counts(1:end),[512 512]);


% Group frequencies into 3x3 squares
x = 2:3:512;
y = 2:3:512;
[xg,yg] = meshgrid(x,y);
subs = reshape(cat(3,xg,yg),[numel(xg) 2]);
inds = sub2ind([512 512],subs(:,1),subs(:,2));
fatmat = padarray(countsIm,[1 1],'symmetric');
averaged = colfilt(fatmat,[3 3],'sliding',@cold);
averaged = averaged(1:512,1:512);
tmp = averaged(inds);
countsIm2 = reshape(tmp,[ceil(512/3) ceil(512/3)]);
expCounts = exp(countsIm2);