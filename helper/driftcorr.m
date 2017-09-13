% Drift Correction

% Use Frame to collect the coordinates of the particles on the frames in
% the interval of interest, want to minimize the for loops
ps = 160; % nm, pixel size
slices = 50;
T = floor(Nframes/slices);
m = 512;
N = zeros(m,m,slices);
% N = zeros(m,m,9);
D = zeros(m,m,slices-1);
driftobjs = [];
for k = 1:size(N,3)
    startframe = (k-1)*T+1;
    endframe = k*T;
    interval = startframe:endframe;
    X = round(objs_link(1,cell2mat(Frame(interval))));
    Y = round(objs_link(2,cell2mat(Frame(interval))));
%     X = round(objs_link(1,Frame{k}),1);
%     Y = round(objs_link(2,Frame{k}),1);
    [N(:,:,k),Xedges,Yedges] = histcounts2(X,Y,m);
    % Perform cross-correlation between each of the projections
    if k == 1
        continue
    else
        tmp = xcorr2(N(:,:,k-1),N(:,:,k));
        D(:,:,k-1) = tmp(257:512+256,257:512+256);
    end
    % Get neighborhoods
    img = D(:,:,k-1);
    tmpobj = analysis_rpstyle(img);
    tmpobj(5,:) = k-1;
    driftobjs = [driftobjs,tmpobj];
end

driftobjs_link = nnlink_rp(driftobjs, 10, 1, true);

% Now everything is linked, get the drift in each direction at each time
% interval
numTrajs = max(driftobjs_link(6,:));
xcoords = driftobjs_link(1,:);
ycoords = driftobjs_link(2,:);
trkID = driftobjs_link(6,:);
XC = NaN(numTrajs,slices); YC = XC;
for k = 1:numTrajs
    currTraj = driftobjs_link(1:2,trkID == k);
    firstframe = min(driftobjs_link(5,trkID == k));
    trajInd = firstframe:firstframe+size(currTraj,2)-1;
    XC(k,trajInd) = currTraj(1,:);
    YC(k,trajInd) = currTraj(2,:);
end

dXC = diff(XC,1,2); dYC = diff(YC,1,2);
xdriftvals = mean(dXC,1,'omitnan')*ps; 
xsem = xdriftvals./sqrt(sum(~isnan(dXC),1));
ydriftvals = mean(dYC,1,'omitnan')*ps;
ysem = ydriftvals./sqrt(sum(~isnan(dYC),1));

frameInt = 2;
timeInt = frameInt*T:frameInt*T:(slices-1)*frameInt*T;
figure,errorbar(timeInt,xdriftvals,xsem,xsem,'b.-')
hold on, errorbar(timeInt,ydriftvals,ysem,ysem,'r.-'), hold off
legend('x','y','Location','best'), xlabel('Time Interval, t (s)'), ylabel('Drift Value (nm)')
title([num2str(slices),' slices']);

% Cubic spline interpolation
timeVec = 2:frameInt:(Nframes-1)*frameInt;
xq = interp1(timeInt,xdriftvals,timeVec,'spline');
yq = interp1(timeInt,ydriftvals,timeVec,'spline');
figure,plot(timeVec,xq,'b.-',timeVec,yq,'r.-')

% Apply correction
% Remove trajectories less than 2 frames
Clong = C(cellfun(@(x) size(x,2),C)>1);
Flong = F(cellfun(@length,F)>1);
% Find location in time
XC = NaN(length(Clong),Nframes); YC = XC;
for k = 1:length(Flong)
    inds = Flong{k};
    XC(k,inds) = Clong{k}(1,:);
    YC(k,inds) = Clong{k}(2,:);
end
% Correct the displacements with the sum of displacement up to that point
xcorr = [NaN(1,T),cumsum(xq(T:end))]; ycorr = [NaN(1,T),cumsum(yq(T:end))];
% xcorr(isnan(xcorr)) = []; ycorr(isnan(ycorr)) = []; 
% XC(isnan(XC)) = []; YC(isnan(YC)) = [];
XCcorr = XC - xcorr/ps; YCcorr = YC - ycorr/ps;
% Plot the cumsum drift correction

Ccorr = cell(size(Clong));
for k = 1:length(Ccorr)
    xs = XCcorr(k,:); ys = YCcorr(k,:); 
    xs(isnan(xs))=[];ys(isnan(ys)) = [];
    % Concatenate the new coordinates with frame numbers
    if ~isempty(xs)
        frames = Flong{k};
        frames = frames(frames>T);
        Ccorr{k} = [xs;ys;frames];
    end  
end

Ccorr = cellfun(@(x) x',Ccorr,'UniformOutput',0);
Ccorrmat = cell2mat(Ccorr);
Clong_col = cellfun(@(x) x',Clong,'UniformOutput',0); 
Flong_col = cellfun(@(x) x',Flong,'UniformOutput',0);
Clongmat = cell2mat(Clong_col); Flongmat = cell2mat(Flong_col);
orig = [Clongmat,Flongmat];

% Display corrected points
figure, h = gca; h.XLim = [1 512]; h.YLim = [1 512]; hold on
for k = T:Nframes
    centers = Ccorrmat(Ccorrmat(:,3) == k,1:2);
    radii = repmat(3,length(centers),1);
    plot(centers(:,1),512-centers(:,2),'g.')
    % viscircles(h,centers,radii,'Color','g');
    origcenters = orig(orig(:,3) == k,1:2);
    radii = repmat(5,length(origcenters),1);
    plot(origcenters(:,1),512-origcenters(:,2),'r.')
    % viscircles(h,origcenters,radii,'Color','r');
    pause(0.05);
end
% Check that the units are right for the original coordinates














