function [Ccorr, Ccorrmat, O] = driftcorr(slices,Nframes,M,objs_link)
% Drift Correction 1 - Correlation Matrix Shift

% Use Frame to collect the coordinates of the particles on the frames in
% the interval of interest, want to minimize the for loops
et = 2; % s, time between frames
T = floor(Nframes/slices);
m = 512; n = 512;
N = zeros(m,m,slices); % Ndil = N;
corrmax = zeros(slices-1,1);

% Keep long trajs only
K = M(cellfun(@length,M)>T);
inds2keep = cat(2,K{:});
objsLinked = objs_link(:,inds2keep);
Frame2 = cell(Nframes,1);
for i = 1:Nframes
    Frame2{i} = find(objsLinked(5,:)==i);
end

Xedges = m/n/2:m/n:n+m/n/2; Yedges = Xedges;
partperint = zeros(slices,1);
for k = 1:slices
    startframe = (k-1)*T+1;
    endframe = k*T;
    interval = startframe:endframe;
    X = round(objsLinked(1,cat(2,Frame2{interval})));
    Y = round(objsLinked(2,cat(2,Frame2{interval}))); 
    N(:,:,k) = histcounts2(X,Y,Xedges,Yedges);
    partperint(k) = sum(sum(N(:,:,k)));
    % Dilate N for better overlap
    % Ndil(:,:,k) = imdilate(N(:,:,k),strel('diamond',1));
    % Perform cross-correlation between each of the projections
    if k == 1
        continue
    else
        tmp = xcorr2(N(:,:,k-1),N(:,:,k));
        D = tmp(m/2+1:m+m/2,m/2+1:m+m/2);
        % Get location of maximum correlation
        % Find first maximum in the image
        corrmax(k-1) = find(D == max(max(D,[],1),[],2),1);
    end
end
% clear tmp D

% Find the difference 
[I,J] = ind2sub([m m],corrmax);
dI = I-m/2; dJ = J-m/2;
Ivect = cumsum(dI); Jvect = cumsum(dJ);
ints = round(linspace(et*T,et*Nframes,slices-1));
timeVec = et:et:Nframes*et;
iq = interp1(ints,Ivect,timeVec,'spline','extrap');
jq = interp1(ints,Jvect,timeVec,'spline','extrap');


% Apply correction to all particles
h = waitbar(0,'Analyzing...');
trkIDall = objs_link(6,:); numAllTrajs = max(trkIDall);
trksAll = cell(numAllTrajs,1); frmsAll = trksAll;
Xall = NaN(numAllTrajs,Nframes); Yall = Xall; bbAll = Xall;
 
for k = 1:numAllTrajs
    trksAll{k} = objs_link(1:2,trkIDall==k);
    frmsAll{k} = objs_link(5,trkIDall==k);
    Xall(k,frmsAll{k}) = objs_link(1,trkIDall==k);
    Yall(k,frmsAll{k}) = objs_link(2,trkIDall==k);
    bbAll(k,frmsAll{k}) = objs_link(3,trkIDall==k); 
    waitbar(k/numAllTrajs);
end
close(h); clear h
XCall = Xall + iq; YCall = Yall + jq;
XCall(XCall<=0) = -1; XCall(XCall>m) = -1;
YCall(YCall<=0) = -1; YCall(YCall>n) = -1;

Ccorr = cell(size(trksAll)); O = Ccorr;
for k = 1:length(Ccorr)
    xs = XCall(k,:); ys = YCall(k,:);
    removedx = xs==-1; removedy = ys==-1;
    todeletex = isnan(xs)|removedx;
    todeletey = isnan(ys)|removedy;
    xs(todeletex|todeletey)=[];ys(todeletex|todeletey) = [];
    if ~isempty(xs)
        frames = frmsAll{k}; 
        remFrames = ismember(frames,find(removedx|removedy));
        molecNum = M{k}; 
        molecNum(remFrames) = [];
        frames(remFrames) = [];
        O{k} = molecNum;
        Ccorr{k} = [xs',ys',frames'];
    end
end
Ccorrmat = cell2mat(Ccorr);
trksmat = cell2mat(cellfun(@(x) x',trksAll,'UniformOutput',0));
frmsmat = cell2mat(cellfun(@(x) x',frmsAll,'UniformOutput',0));
orig = [trksmat,frmsmat];

% Display corrected points
figure, h = gca; h.XLim = [1 512]; h.YLim = [1 512]; hold on
for k = 1:Nframes
    origcenters = orig(orig(:,3) == k,1:2);
    plot(origcenters(:,1),512-origcenters(:,2),'r.')
    centers = Ccorrmat(Ccorrmat(:,3) == k,1:2);
    plot(centers(:,1),512-centers(:,2),'b.')
    pause(0.025)
end
% figure, h = gca; h.XLim = [1 512]; h.YLim = [1 512]; hold on
% plot(objs_link(1,:),512-objs_link(2,:),'r.')
% plot(XCall,512-YCall,'b.'), hold on