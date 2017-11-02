% Load image
iminf = imfinfo(filename);
N = length(iminf);
m = iminf(1).Width; n = iminf(1).Height;
img = zeros(m,n,N,'uint16');
TifLink = Tiff(filename,'r');
for k = 1:N
TifLink.setDirectory(k);
img(:,:,k) = TifLink.read();
end
TifLink.close();

% Get standard deviation projection
stdim = std(double(img),[],3);

% Convert image to double range
I1 = stdim/max(stdim(:));

% Binarize with Otsu's method to generate seeding locations 
I2 = imbinarize(I1,'adaptive');
start = regionprops('table',I2,'Centroid');
start = start.Centroid;

% % Find mask overlaps
% traj = cell(length(CC.PixelIdxList),1); tol = 0.71;
% for k = 1:length(CC.PixelIdxList)
%     [I,J] = ind2sub([512 512],CC.PixelIdxList{k});
%     inds = find(objs_link(1,:) < max(I)+tol & objs_link(1,:) > min(I)-tol...
%         & objs_link(2,:) < max(J)+tol & objs_link(2,:) > min(J)-tol);
%     if ~isempty(inds)
%         traj{k} = [inds; objs_link([1,2,3,5,6],inds)];
%     end
% end

% Find the objects within tolerated distance of the seeded points using
% kmeans clustering
tol = 3;
coords = [Ccorrmat(:,1), Ccorrmat(:,2)];
[idx,c,sumd,d] = kmeans(coords,[],'Start',start);
tocluster = d; tocluster(d>tol) = NaN;

% Only allow one cluster assignment for each object. Take minimum
% distance
numclose = sum(~isnan(tocluster),2);
if sum(numclose > 1) > 0
    for k = find(sum(numclose > 1))
        mindist = tocluster(k,:) == min(tocluster(k,:));
        tocluster(k,~mindist) = NaN;
    end
end

% Organize the clusters into object trajectories. Objtrajtime is separated
% in time by NaNs for tracing purposes.
trajnums = cell(length(Ccorr),1);
for k = 1:length(Ccorr)
    trajnums{k} = k*ones(size(Ccorr{k},1),1);
end
trajnumsmat = cell2mat(trajnums);

objtrajinf = cell(size(tocluster,2)); objtrajtime = objtrajinf; 
spotEvents = struct;
for k = 1:size(tocluster,2)
    trajtoinc = trajnummat(~isnan(tocluster(:,k)));
    objnum = [R(trajtoinc).trajectory];
    tmpinf(1:2,:) = [Ccorrmat(~isnan(tocluster(:,k)))'; objs_link([3;5;6],objnum)];
    [frmsrtd,sortind] = sort(tmpinf(3,:),2,'ascend');
    objtrajinf{k} = tmpinf(:,sortind);
    tmptraj = NaN(1,N);
    tmptraj(frmsrted) = tmpinf(:,sortind);
    objtrajtime{k} = tmptraj;
    spotEvents(k).trajectory = objnum(sortind);
end