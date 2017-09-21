% Drift Correction

% Use Frame to collect the coordinates of the particles on the frames in
% the interval of interest, want to minimize the for loops
ps = 160; % nm, pixel size
et = 2; % s, exposure time
slices = 10;
T = floor(Nframes/slices);
m = 512;
N = zeros(m,m,slices); Ndil = N;
% N = zeros(m,m,9);
D = zeros(m,m,slices-1);
driftobjs = [];
corrmax = zeros(slices-1,1);
% try long trajs only
K = M(cellfun(@length,M)>T);
inds2keep = cat(2,K{:});
objsLinked = objs_link(:,inds2keep);
Frame2 = cell(Nframes,1);
for i = 1:Nframes
    Frame2{i} = find(objsLinked(5,:)==i);
end

Xedges = 0.5:1:m+0.5; Yedges = Xedges;
partperint = zeros(slices,1);
for k = 1:slices
    startframe = (k-1)*T+1;
    endframe = k*T;
    interval = startframe:endframe;
    X = round(objsLinked(1,cat(2,Frame2{interval})));
    Y = round(objsLinked(2,cat(2,Frame2{interval}))); % includes m-numbers not in objsLinked
    N(:,:,k) = histcounts2(X,Y,Xedges,Yedges);
    partperint(k) = sum(sum(N(:,:,k)));
    % Dilate N for better overlap
    % Ndil(:,:,k) = imdilate(N(:,:,k),strel('diamond',1));
    % Perform cross-correlation between each of the projections
    if k == 1
        continue
    else
        tmp = xcorr2(N(:,:,k-1),N(:,:,k));
        D(:,:,k-1) = tmp(m/2+1:m+m/2,m/2+1:m+m/2);
        % Get location of maximum correlation
        corrmax(k-1) = find(D(:,:,k-1) == max(max(D(:,:,k-1),[],1),[],2));
    end
    
    % Get neighborhoods
%     img = D(:,:,k-1);
%     tmpobj = analysis_rpstyle(img);
%     tmpobj(5,:) = k-1;
%     driftobjs = [driftobjs,tmpobj];
end

% Find the difference 
[I,J] = ind2sub([m m],corrmax);
dI = I-m/2; dJ = J-m/2;
Ivect = cumsum(dI); Jvect = cumsum(dJ);
ints = round(linspace(et*T,et*Nframes,slices-1));
% figure,plot(ints,dI,'k.-',ints,dJ,'b^-')
figure,plot(ints,Ivect,'k.-',ints,Jvect,'b^-')
title([num2str(slices),' slices, ',num2str(mean(partperint)),' particles per interval'])
timeVec = et:et:Nframes*et;
iq = interp1(ints,Ivect,timeVec,'spline');
jq = interp1(ints,Jvect,timeVec,'spline');
figure,plot(timeVec,iq,'k.-',timeVec,jq,'b^-')

% Apply correction
trkID = objsLinked(6,:); numTrajs = max(trkID); 
trks = cell(numTrajs,1); frms = trks;
XC = NaN(length(trks),Nframes); YC = XC;
for k = 1:numTrajs
    trks{k} = objsLinked(1:2,trkID==k);
    frms{k} = objsLinked(5,trkID==k);
    XC(k,frms{k}) = objsLinked(1,trkID==k);
    YC(k,frms{k}) = objsLinked(2,trkID==k);
end
XCcorr = XC + iq.*ps; YCcorr = YC + jq.*ps;

% Convert into individual trajectory cell form
Ccorr = cell(size(trks));
for k = 1:length(Ccorr)
    xs = XCcorr(k,:); ys = YCcorr(k,:);
    xs(isnan(xs))=[];ys(isnan(ys)) = [];
    if ~isempty(xs)
        frames = frms{k};
        Ccorr{k} = [xs',ys',frames'];
    end
end
Ccorrmat = cell2mat(Ccorr);
trksmat = cell2mat(cellfun(@(x) x',trks,'UniformOutput',0));
frmsmat = cell2mat(cellfun(@(x) x',frms,'UniformOutput',0));
orig = [trksmat,frmsmat];

% Display corrected points
figure, h = gca; h.XLim = [1 512]; h.YLim = [1 512]; hold on
for k = 1:Nframes
    centers = Ccorrmat(Ccorrmat(:,3) == k,1:2);
    plot(centers(:,1),512-centers(:,2),'g.')
    % radii = repmat(3,length(centers),1);
    % viscircles(h,centers,radii,'Color','g');
    origcenters = orig(orig(:,3) == k,1:2);
    plot(origcenters(:,1),512-origcenters(:,2),'r.')
    % radii = repmat(5,length(origcenters),1);
    % viscircles(h,origcenters,radii,'Color','r');
    pause(0.05);
end

% Ccorr = cell(size(Clong));
% for k = 1:length(Ccorr)
%     xs = XCcorr(k,:); ys = YCcorr(k,:); 
%     xs(isnan(xs))=[];ys(isnan(ys)) = [];
%     % Concatenate the new coordinates with frame numbers
%     if ~isempty(xs)
%         frames = Flong{k};
%         frames = frames(frames>T);
%         Ccorr{k} = [xs;ys;frames];
%     end  
% end
% 
% Ccorr = cellfun(@(x) x',Ccorr,'UniformOutput',0);
% Ccorrmat = cell2mat(Ccorr);
% Clong_col = cellfun(@(x) x',Clong,'UniformOutput',0); 
% Flong_col = cellfun(@(x) x',Flong,'UniformOutput',0);
% Clongmat = cell2mat(Clong_col); Flongmat = cell2mat(Flong_col);
% orig = [Clongmat,Flongmat];
    
    
% 
% 
% driftobjs_link = nnlink_rp(driftobjs, 10, 1, true);
% 
% % Now everything is linked, get the drift in each direction at each time
% % interval
% numTrajs = max(driftobjs_link(6,:));
% xcoords = driftobjs_link(1,:);
% ycoords = driftobjs_link(2,:);
% trkID = driftobjs_link(6,:);
% XC = NaN(numTrajs,slices); YC = XC;
% for k = 1:numTrajs
%     currTraj = driftobjs_link(1:2,trkID == k);
%     firstframe = min(driftobjs_link(5,trkID == k));
%     trajInd = firstframe:firstframe+size(currTraj,2)-1;
%     XC(k,trajInd) = currTraj(1,:);
%     YC(k,trajInd) = currTraj(2,:);
% end
% 
% dXC = diff(XC,1,2); dYC = diff(YC,1,2);
% xdriftvals = mean(dXC,1,'omitnan')*ps; 
% xsem = xdriftvals./sqrt(sum(~isnan(dXC),1));
% ydriftvals = mean(dYC,1,'omitnan')*ps;
% ysem = ydriftvals./sqrt(sum(~isnan(dYC),1));
% 
% frameInt = 2;
% timeInt = frameInt*T:frameInt*T:(slices-1)*frameInt*T;
% figure,errorbar(timeInt,xdriftvals,xsem,xsem,'b.-')
% hold on, errorbar(timeInt,ydriftvals,ysem,ysem,'r.-'), hold off
% legend('x','y','Location','best'), xlabel('Time Interval, t (s)'), ylabel('Drift Value (nm)')
% title([num2str(slices),' slices']);
% 
% % Cubic spline interpolation
% timeVec = 2:frameInt:(Nframes-1)*frameInt;
% xq = interp1(timeInt,xdriftvals,timeVec,'spline');
% yq = interp1(timeInt,ydriftvals,timeVec,'spline');
% figure,plot(timeVec,xq,'b.-',timeVec,yq,'r.-')
% 
% % Apply correction
% % Remove trajectories less than 2 frames
% Clong = C(cellfun(@(x) size(x,2),C)>1);
% Flong = F(cellfun(@length,F)>1);
% % Find location in time
% XC = NaN(length(Clong),Nframes); YC = XC;
% for k = 1:length(Flong)
%     inds = Flong{k};
%     XC(k,inds) = Clong{k}(1,:);
%     YC(k,inds) = Clong{k}(2,:);
% end
% % Correct the displacements with the sum of displacement up to that point
% xcorr = [NaN(1,T),cumsum(xq(T:end))]; ycorr = [NaN(1,T),cumsum(yq(T:end))];
% % xcorr(isnan(xcorr)) = []; ycorr(isnan(ycorr)) = []; 
% % XC(isnan(XC)) = []; YC(isnan(YC)) = [];
% XCcorr = XC - xcorr/ps; YCcorr = YC - ycorr/ps;
% % Plot the cumsum drift correction
% 
% Ccorr = cell(size(Clong));
% for k = 1:length(Ccorr)
%     xs = XCcorr(k,:); ys = YCcorr(k,:); 
%     xs(isnan(xs))=[];ys(isnan(ys)) = [];
%     % Concatenate the new coordinates with frame numbers
%     if ~isempty(xs)
%         frames = Flong{k};
%         frames = frames(frames>T);
%         Ccorr{k} = [xs;ys;frames];
%     end  
% end
% 
% Ccorr = cellfun(@(x) x',Ccorr,'UniformOutput',0);
% Ccorrmat = cell2mat(Ccorr);
% Clong_col = cellfun(@(x) x',Clong,'UniformOutput',0); 
% Flong_col = cellfun(@(x) x',Flong,'UniformOutput',0);
% Clongmat = cell2mat(Clong_col); Flongmat = cell2mat(Flong_col);
% orig = [Clongmat,Flongmat];
% 
% % Display corrected points
% figure, h = gca; h.XLim = [1 512]; h.YLim = [1 512]; hold on
% for k = T:Nframes
%     centers = Ccorrmat(Ccorrmat(:,3) == k,1:2);
%     radii = repmat(3,length(centers),1);
%     plot(centers(:,1),512-centers(:,2),'g.')
%     % viscircles(h,centers,radii,'Color','g');
%     origcenters = orig(orig(:,3) == k,1:2);
%     radii = repmat(5,length(origcenters),1);
%     plot(origcenters(:,1),512-origcenters(:,2),'r.')
%     % viscircles(h,origcenters,radii,'Color','r');
%     pause(0.05);
% end
% % Check that the units are right for the original coordinates
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
