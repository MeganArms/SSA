% % Load image
% iminf = imfinfo(filename);
% N = length(iminf);
% m = iminf(1).Width; n = iminf(1).Height;
% img = zeros(m,n,N,'uint16');
% TifLink = Tiff(filename,'r');
% for k = 1:N
% TifLink.setDirectory(k);
% img(:,:,k) = TifLink.read();
% end
% TifLink.close();
% 
% % Get standard deviation projection
% stdim = std(double(img),[],3);
% 
% % Convert image to double range
% I1 = stdim/max(stdim(:));
% 
% % Binarize with Otsu's method to generate seeding locations 
% I2 = imbinarize(I1,'adaptive');
% start = regionprops('table',I2,'Centroid');
% start = start.Centroid;

% Find the objects within tolerated distance of the seeded points using
% kd-tree algorithm
tol = 3;
idx = rangesearch(coords,startpad,tol);

% Organize the clusters into object trajectories. spotEvents is separated
% in time by NaNs for tracing purposes.
trajnums = cell(length(Ccorr),1);
for k = 1:length(Ccorr)
    trajnums{k} = k*ones(size(Ccorr{k},1),1);
end
trajnumsmat = cell2mat(trajnums);

loc = ~cellfun(@isempty,idx);
inds = find(loc);
objtrajinf = cell(sum(loc),1);  
spotEvents = struct; 
for k = 1:length(inds)
    trajtoinc = unique(trajnumsmat(idx{k}));
    l = inds(k);
    if ~isempty(trajtoinc)
        objnum = [R(trajtoinc).trajectory];
        tmpcoords = cell2mat(Ccorr(trajtoinc));
        tmpinf = [tmpcoords(:,1:2)'; objs_link([3;5;6],objnum)];
        [frmsrtd,sortind] = sort(tmpinf(4,:),2,'ascend');
        objtrajinf{k} = tmpinf(:,sortind);
        tmptraj = NaN(1,N); tmptraj(frmsrtd) = objnum(sortind);
        spotEvents(k).trajectory = tmptraj;
    end
end