% % Drift Correction
% 
% % Use Frame to collect the coordinates of the particles on the frames in
% % the interval of interest, want to minimize the for loops
% T = floor(Nframes/10);
% m = 512;
% N = zeros(m,m,10);
% % N = zeros(m,m,9);
% D = zeros(m,m,size(N,3)-1);
% driftobjs = [];
% for k = 1:size(N,3)
%     startframe = (k-1)*T+1;
%     endframe = k*T;
%     interval = startframe:endframe;
%     X = round(objs_link(1,cell2mat(Frame(interval))),1);
%     Y = round(objs_link(2,cell2mat(Frame(interval))),1);
% %     X = round(objs_link(1,Frame{k}),1);
% %     Y = round(objs_link(2,Frame{k}),1);
% %     [N(:,:,k),Xedges,Yedges] = histcounts2(X,Y,m);
%     % Perform cross-correlation between each of the projections
%     if k ~= 1
%         tmp = xcorr2(N(:,:,k-1),N(:,:,k));
%         D(:,:,k-1) = tmp(257:512+256,257:512+256);
%     else
%         continue
%     end
%     % Get neighborhoods
%     img = D(:,:,k-1);
%     driftobjs = [driftobjs,analysis_rpstyle(img)];
% end
% 
% driftobj_link = nnlink_rp(driftobjs, 10, 1, true);

% Now everything is linked, get the drift in each direction at each time
% interval
numTrajs = max(objs_link(6,:));
xcoords = objs_link(1,:);
ycoords = objs_link(2,:);
trkID = objs_link(6,:);
XC = cell(numTrajs,1); YC = XC;
for k = 1:numTrajs
    XC{k} = objs_link(1,trkID == k);
    YC{k} = objs_link(2,trkID == k);
end

dXC = cellfun(@diff,XC,'UniformOutput',0); % This isn't right! is it?
% yeah it is because I want the drift between points, not the points
% themselves
dYC = cellfun(@diff,YC,'UniformOutput',0);
lengths = cellfun(@length,XC);
longest = max(lengths);

dXCpadded = cellfun(@(x) padarray(x,[0 longest-length(x)],NaN,'post'),dXC,'UniformOutput',0);
dYCpadded = cellfun(@(x) padarray(x,[0 longest-length(x)],NaN,'post'),dYC,'UniformOutput',0);
dXCmat = cat(1,dXCpadded{:}); YCmat = cat(1,dYCpadded{:});

xdriftvals = mean(dXCmat,1,'omitnan')*160;
ydriftvals = mean(YCmat,1,'omitnan')*160;

timeInt = 200:200:longest*200;
figure,plot(timeInt,xdriftvals,'b.-'),hold on, plot(timeInt,ydriftvals,'r.-'),hold off
legend('x','y','Location','best'), xlabel('Time Interval, t (s)'), ylabel('Drift Value (nm)')

% Cubic spline interpolation

