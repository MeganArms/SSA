% Drift Correction 2 - center of mass shift

% Get molecular indices of particles on the frame
dCOM = zeros(Nframes-1,2);
% h = waitbar(0,'Analyzing...');
% figure, h = gca; h.XLim = [1 512]; h.YLim = [1 512]; hold on
for k = 2:Nframes
    traj2use = cellfun(@(x) ismember(k-1,x),F) & ...
        cellfun(@(x) ismember(k,x),F);
    inds = 1:length(C); indskept = inds(traj2use);
    traj2use = C(traj2use); 
    coords2use1 = cell(length(traj2use)); coords2use2 = coords2use1;
    for l = 1:length(traj2use)
        p = indskept(l);
        coords2use1{l} = traj2use{l}(:,F{p}==k-1);
        coords2use2{l} = traj2use{l}(:,F{p}==k);
    end
    coords2use1 = cat(2,coords2use1{:}); coords2use2 = cat(2,coords2use2{:});
    COM1 = mean(coords2use1,2,'omitnan'); COM2 = mean(coords2use2,2,'omitnan');
    dCOM(k-1,:) = [COM2(1)-COM1(1),COM2(2)-COM1(2)];
%     plot(coords2use1(1,:),512-coords2use1(2,:),'k.'), 
%     plot(coords2use2(1,:),512-coords2use2(2,:),'b.'),
%     plot(COM(k,1),512-COM(k,2),'r.'), 
%     pause(0.025);
    waitbar(k/Nframes);
end
cumDiff = [[0;0],cumsum(dCOM,1)'];

% close(h)

% Apply correction
trkID = objs_link(6,:); numTrajs = max(trkID); 
trks = cell(numTrajs,1); frms = trks;
XC = NaN(length(trks),Nframes); YC = XC;
for k = 1:numTrajs
    trks{k} = objs_link(1:2,trkID==k);
    frms{k} = objs_link(5,trkID==k);
    XC(k,frms{k}) = objs_link(1,trkID==k);
    YC(k,frms{k}) = objs_link(2,trkID==k);
end
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
