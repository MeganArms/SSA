% Drift Correction 1 - Correlation Matrix Shift

% Use Frame to collect the coordinates of the particles on the frames in
% the interval of interest, want to minimize the for loops
ps = 160; % nm, pixel size
et = 2; % s, exposure time
slices = 10;
T = floor(Nframes/slices);
m = 512; n = 512;
N = zeros(m,m,slices); % Ndil = N;
% D = zeros(m,m,slices-1);
corrmax = zeros(slices-1,1);
I = zeros(slices-1,1); J = zeros(slices-1,1);
Ilo = I; Iup = I; Jlo = J; Jup = J;

% Keep long trajs only
K = M(cellfun(@length,M)>T);
inds2keep = cat(2,K{:});
objsLinked = objs_link(:,inds2keep);
Frame2 = cell(Nframes,1);
for i = 1:Nframes
    Frame2{i} = find(objsLinked(5,:)==i);
end

h = waitbar(0,'Analyzing...');
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
    waitbar(k/slices);
end
close(h); clear h

% Find the difference 
[I,J] = ind2sub([m m],corrmax);
dI = I-m/2; dJ = J-m/2;
Ivect = cumsum(dI); Jvect = cumsum(dJ);
% Ivect = cumsum(I); Jvect = cumsum(J);
% Iupvect = (Iup-m/2 - dI) + Ivect; Ilovect = Ivect - (dI - (Ilo-m/2));
% Jupvect = (Jup-m/2 - dJ) + Jvect; Jlovect = Jvect - (dJ - (Jlo-m/2)); 
ints = round(linspace(et*T,et*Nframes,slices-1));
% figure,plot(ints,I,'k.-',ints,J,'b^-')
% figure,plot(ints,Ivect,'k.-',ints,Iupvect,'k.--',ints,Ilovect,'k.--',...
%     ints,Jvect,'b^-',ints,Jupvect,'b^--',ints,Jlovect,'b^--')
% title([num2str(slices),' slices, ',num2str(mean(partperint)),' particles per interval'])
timeVec = et:et:Nframes*et;
% iloq = interp1(ints,Ilovect,timeVec,'linear','extrap');
iq = interp1(ints,Ivect,timeVec,'spline','extrap');
% iupq = interp1(ints,Iupvect,timeVec,'linear','extrap');
% jloq = interp1(ints,Jlovect,timeVec,'linear','extrap');
jq = interp1(ints,Jvect,timeVec,'spline','extrap');
% jupq = interp1(ints,Jupvect,timeVec,'linear','extrap');
% figure,plot(timeVec,iq,'k.-',timeVec,iloq,'k--',timeVec,iupq,'k--',...
%     timeVec,jq,'b^-',timeVec,jloq,'b--',timeVec,jupq,'b--')
% figure, plot(timeVec,iq,'k.-',timeVec,jq,'b^-')

% Apply correction
trkID = objsLinked(6,:); numTrajs = max(trkID); 
trks = cell(numTrajs,1); frms = trks;
XC = NaN(length(trks),Nframes); YC = XC; bb = XC;
for k = 1:numTrajs
    trks{k} = objsLinked(1:2,trkID==k);
    frms{k} = objsLinked(5,trkID==k);
    XC(k,frms{k}) = objsLinked(1,trkID==k);
    YC(k,frms{k}) = objsLinked(2,trkID==k);
    bb(k,frms{k}) = objsLinked(3,trkID==k);
end
XCcorr = XC + iq; YCcorr = YC + jq;

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
    origcenters = orig(orig(:,3) == k,1:2);
    plot(origcenters(:,1),512-origcenters(:,2),'r.')
    centers = Ccorrmat(Ccorrmat(:,3) == k,1:2);
    plot(centers(:,1),512-centers(:,2),'b.')
    pause(0.025)
end
