function [spotEvents3, eventfreq] = eventlinks9(tol,Ccorr,O,B,varargin)

%% Find locations of interest
if ~isempty(varargin)
    F = varargin{1};
    Ccorr = cellfun(@(x)x',Ccorr,'UniformOutput',false);
end

longCoords = Ccorr; clengths = zeros(length(Ccorr),1);
objnum = O;
seedCoords = zeros(length(Ccorr),2);
for k = 1:length(longCoords)
    % clengths(k) = size(longCoords{k},2);
    clengths(k) = size(longCoords{k},1);
    if clengths(k) <= 1
        longCoords{k} = [];
        objnum{k} = [];
    else
        seedCoords(k,:) = mean(longCoords{k}(:,1:2),1);
        % seedCoords(k,:) = mean(longCoords{k},2)';
    end
end

inds = 1:length(Ccorr);
inds = inds(clengths>1);
seedCoords = seedCoords(clengths>1,:);

eucMatrix = triu(squareform(pdist(seedCoords)));
nn = eucMatrix < tol & eucMatrix > 0;

%% Find trajectories in the same location
combineIdx = cell(length(seedCoords),1);
beenadded = [];
for k = 1:length(seedCoords)
    toadd = inds(nn(k,:));
    if sum(ismember([inds(k),toadd],beenadded))==0
        combineIdx{k} = [inds(k),toadd];
        nn(nn(k,:),:) = 0;
        beenadded = [beenadded,combineIdx{k}];
    end
end

%% Create new, connected trajectories
spotEvents = struct; 
eventfreq = zeros(length(combineIdx),1);
l = 1;
for k = 1:length(combineIdx)
    objnums = [objnum{combineIdx{k}}];
    brightness = [B{combineIdx{k}}];
    frms = [F{combineIdx{k}}];
    if ~isempty(objnums)
        tocmbn = cat(1,longCoords{combineIdx{k}});
        
        % Check for bad connections from RP
        coord2check = tocmbn(:,1:2);
        disttest = triu(squareform(pdist(coord2check)));
        bb = disttest(1,:) < 4*tol;
%         A = diag(diag(disttest,1),1);
%         [bbi,bbj] = ind2sub(size(disttest),find(A < tol & A > 0));
%         bb = unique([bbi;bbj]);
        objnums = objnums(bb);
        tocmbn = tocmbn(bb,:);
        brightness = brightness(bb);
        
        % Make new connections
%         frames = round(tocmbn(:,3));
        frames = frms(bb);
        tmptraj = NaN(1,max(frames)-min(frames)+1);
        [frmsrtd,srtind] = sort(frames,'ascend');
        indnum = frmsrtd-min(frmsrtd)+1;
        tmptraj(indnum) = objnums(srtind);
        spotEvents(l).trajectory = tmptraj;
        spotEvents(l).coordinates = tocmbn(srtind,1:2);
        spotEvents(l).brightness = brightness(srtind);
        spotEvents(l).frames = frmsrtd;
        CC = bwconncomp(isnan(tmptraj));
        lengths = cellfun(@length,CC.PixelIdxList);
        eventfreq(l) = sum(lengths>1)+1; % Dark time must be greater than one frame
        spotEvents(l).eventfreq = eventfreq(l);
        spotEvents(l).std = sqrt(sum(std(tocmbn(:,1:2)).^2));
        l = l + 1;
    end
end
eventfreq(eventfreq==0) = [];

%% Cleanup

spotEvents3 = struct; l = 1;
for k = 1:length(spotEvents)
    if spotEvents(k).std < 4
        spotEvents3(l).trajectory = spotEvents(k).trajectory;
        spotEvents3(l).coordinates = spotEvents(k).coordinates;
        spotEvents3(l).brightness = spotEvents(k).brightness;
        spotEvents3(l).frames = spotEvents(k).frames;
        spotEvents3(l).eventfreq = spotEvents(k).eventfreq;
        spotEvents3(l).std = spotEvents(k).std;
        l = l + 1;
    end
end
