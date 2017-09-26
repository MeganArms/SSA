% Drift Correction 2 - center of mass shift

% Get molecular indices of particles on the frame
T = Nframes;
inds = 1:length(C); indskept = inds(cellfun(@length,C)>=T);
c = C(indskept);
K = M(indskept);
inds2keep = cat(2,K{:});
objsLinked = objs_link(:,inds2keep);
COM = zeros(Nframes,2);
% h = waitbar(0,'Analyzing...');
% figure, h = gca; h.XLim = [1 512]; h.YLim = [1 512]; hold on
for k = 1:Nframes
%     if k == Nframes
%         traj2use = cellfun(@(x) ismember(k-1,x),F) & ...
%             cellfun(@(x) ismember(k,x),F);
%     else
%         traj2use = cellfun(@(x) ismember(k,x),F) & ...
%             cellfun(@(x) ismember(k+1,x),F);
%     end
%     inds = 1:length(C); indskept = inds(traj2use);
%     traj2use = C(traj2use); 
%     coords2use = cell(length(traj2use));
%     for l = 1:length(traj2use)
%         p = indskept(l);
%         coords2use{l} = traj2use{l}(:,F{p}==k);
%     end
    coords2use = cell(length(c),1);
    lengths = zeros(length(c),1);
    for l = 1:length(c)
        p = indskept(l);
        coords2use{l} = c{l}(:,F{p}==k);
        if sum(F{p}==k)
            lengths(l) = length(c{l});
        end
    end
    lengths = lengths(lengths~=0);
    % coords2use = cellfun(@(x) x(:,k),traj2use,'UniformOutput',0);
    coords2use = cat(2,coords2use{:});
    COM(k,:) = [mean(coords2use(1,:)),mean(coords2use(2,:))];
    % COM(k,:) = [sum(coords2use(1,:)./lengths'),sum(coords2use(2,:)./lengths')];
%     plot(coords2use(1,:),512-coords2use(2,:),'k.'), 
%     plot(COM(k,1),512-COM(k,2),'r.'), 
%     pause(0.025);
    % waitbar(k/Nframes);
end
dCOM = diff(COM,1,1);
cumDiff = [[0;0],cumsum(dCOM,1)'];

% close(h)
% figure
% plot(timeVec(1:end-1),dCOM(:,1)','k^-',timeVec(1:end-1),dCOM(:,2)','b.-');


% Apply correction
trkID = objsLinked(6,:); numTrajs = length(c); 
trks = cell(numTrajs,1); frms = trks;
XC = NaN(length(trks),Nframes); YC = XC;
for k = 1:numTrajs
    trks{k} = objsLinked(1:2,trkID==k);
    frms{k} = objsLinked(5,trkID==k);
    XC(k,frms{k}) = objsLinked(1,trkID==k);
    YC(k,frms{k}) = objsLinked(2,trkID==k);
end
XCcorr = XC - cumDiff(1,:); YCcorr = YC - cumDiff(2,:);

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
    pause(0.025);
end
