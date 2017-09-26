% Drift Correction 2 - center of mass shift

% Convert cell indexing to matrices of coordinates
trkID = objs_link(6,:); numTrajs = max(trkID); 
trks = cell(numTrajs,1); frms = trks;
XC = NaN(length(trks),Nframes); YC = XC;
for k = 1:numTrajs
    trks{k} = objs_link(1:2,trkID==k);
    frms{k} = objs_link(5,trkID==k);
    XC(k,frms{k}) = objs_link(1,trkID==k);
    YC(k,frms{k}) = objs_link(2,trkID==k);
end

% Get center of mass displacements
% h = waitbar(0,'Analyzing...');
% dx = zeros(numTrajs,Nframes-1); dy = dx;
dCOM = zeros(Nframes-1,2);
for k = 2:Nframes
    ind = ~isnan(XC(:,k-1)) & ~isnan(XC(:,k));
%     dx(k-1,:) = diff([XC(ind,k-1),XC(ind,k)],1,2);
%     dy(k-1,:) = diff([YC(ind,k-1),YC(ind,k)],1,2);
    COM1 = [mean(XC(ind,k-1)),mean(YC(ind,k-1))];
    COM2 = [mean(XC(ind,k)),mean(YC(ind,k))];
    dCOM(k-1,:) = [COM2(1)-COM1(1),COM2(2)-COM1(2)];
    waitbar(k/Nframes)
end
cumDiff = [[0;0],cumsum(dCOM,1)'];
% close(h)

% Apply correction
XCcorr = XC - cumDiff(1,:); YCcorr = YC - cumDiff(2,:);

% Convert into individual trajectory cell form
Ccorr = cell(size(trks));
for k = 1:length(Ccorr)
    xs = XCcorr(k,:); ys = YCcorr(k,:);
    xs(isnan(xs)) = []; ys(isnan(ys)) = [];
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
    origcenters = orig(orig(:,3) == k,1:2);
    plot(origcenters(:,1),512-origcenters(:,2),'r.')
    pause(0.025);
end
