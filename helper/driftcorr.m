% Drift Correction

% Use Frame to collect the coordinates of the particles on the frames in
% the interval of interest, want to minimize the for loops
ps = 160; % nm, pixel size
slices = 10;
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
% Correct the displacements
dXC = diff(XC,1,2)*ps; dYC = diff(YC,1,2)*ps;
xcorr = [NaN(1,T-1),cumsum(xq(T:end))]; ycorr = [NaN(1,T-1),cumsum(yq(T:end))];
XCcorr = XC - xcorr; YCcorr = YC - ycorr;
% dXCcorr = dXC - xq; dYCcorr = dYC - yq;
% Convert displacements back into coordinates - this only corrects drift
% after the initial point. Need to correct it over the whole thing, i.e.
% even the first point in a trajectory that appears in the middle of the
% video has experienced drift.
firstind = find(~isnan(XC(k,:)),1);
firstcoord = XC(k,firstind);
XCcorr(k,firstind:end) = [firstcoord, firstcoord + cumsum(dXCcorr(k,firstind:end,'omitnan'))];

XCcorr = XC + xcorr;
















