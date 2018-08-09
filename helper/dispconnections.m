function dispconnections(I,ctrajs)

% Display the connections made between dark trajectories.

if ischar(I)
    % Load image
    iminf = imfinfo(filename);
    N = length(iminf);
    m = iminf(1).Width; n = iminf(1).Height;
    img = zeros(m,n,N,'uint16');
    TifLink = Tiff(filename,'r');
    for k = 1:N
        TifLink.setDirectory(k);
        img(:,:,k) = TifLink.read();
        w = warning('query','last');
        id = w.identifier;
        warning('off',id);
    end
    TifLink.close();
else
    N = size(I,3);
end

numframes = zeros(length(ctrajs),1);
for k = 1:length(ctrajs)
    numframes(k) = ctrajs(k).frames(end)-ctrajs(k).frames(1)+1;
end
% Generate map of connection events
% avgposns = cell(length(ctrajs),1);
sepTrajs = cell(length(ctrajs),1);
avgposc = zeros(sum(numframes),6);
p = 1;
for k = 1:length(ctrajs)
    isnans = isnan(ctrajs(k).trajectory);
    CC = bwconncomp(isnans);
    pixidxlist = CC.PixelIdxList;
    idxlist = pixidxlist(cellfun(@length,pixidxlist)>1);
    inds = cellfun(@(x)[x(1)-1,x(end)+1],idxlist,'UniformOutput',false);
    inds = [1,[inds{:}]];
    coords = nan(length(isnans),2);
    frames = nan(length(isnans),1);
    frames(ctrajs(k).frames-ctrajs(k).frames(1)+1) = ctrajs(k).frames;
    coords(ctrajs(k).frames-ctrajs(k).frames(1)+1,1) = ctrajs(k).coordinates(:,1);
    coords(ctrajs(k).frames-ctrajs(k).frames(1)+1,2) = ctrajs(k).coordinates(:,2);
    frames2 = frames;
    frames(isnan(frames)) = frames([isnan(frames(2:end));false])+1;
    coords1 = coords([isnan(frames2(2:end));false],:);
    coords2 = coords([false;isnan(frames2(1:end-1))],:);
    coordstack = cat(3,coords1,coords2);
    replcoords = mean(coordstack,3);
    coords(isnan(frames2),:) = replcoords;
    cmp = colormap('parula');
    reps = ceil(length(ctrajs)/length(cmp));
    cmap = repmat(cmp,reps,1);
    cmap = cmap(1:length(ctrajs),:);
    % avgposn1 = zeros(length(inds)/2,2);
    % Add frame # to skipped frames and avg coords to skipped frames
    minitraj = zeros(numframes(k),6);
    q = 1;
    for l = 1:2:length(inds)
        if l == length(inds)
            indice = inds(end);
            minitraj(q:q+numframes(k)-inds(end),:) = [coords(indice:end,:),frames(indice:end),repmat(cmap(k,:),length(coords)-indice+1,1)];
        else
            indices = inds(l:l+1);
            % avgposn1(l,:) = mean(coords(indices,:),1);
            minitraj(q:q+diff(indices),:) = [coords(indices(1):indices(2),:),frames(indices(1):indices(2)),repmat(cmap(k,:),diff(indices)+1,1)];
            q = q + diff(indices) + 1;
        end
    end
    % avgposns{k} = avgposn1;
    sepTrajs{k} = minitraj;
    avgcoords = mean(coords,1,'omitnan');
    avgposc(p:p+length(ctrajs(k).frames)-1,:) = [repmat(avgcoords,length(ctrajs(k).frames),1),ctrajs(k).frames',repmat([1 0 0],length(ctrajs(k).frames),1)];
    p = p + length(ctrajs(k).frames);
end
disconnectedtrajs = cell2mat(sepTrajs);
% Display the trajectories
centers = nan(length(disconnectedtrajs),2);
colors = nan(length(disconnectedtrajs),3);
avgcenters = nan(length(avgposc),2);
avgcolors = nan(length(avgposc),3);
figure
for k = 1:N
    imshow(imadjust(I(:,:,k))), hold on
    title('Each connected traj is a separate color')
    centers1 = disconnectedtrajs(disconnectedtrajs(:,3)==k,1:2);
    colors1 = disconnectedtrajs(disconnectedtrajs(:,3)==k,4:6);
    centers(l:l+length(centers1)-1,:) = centers1;
    colors(l:l+length(colors1)-1,:) = colors1;
    l = l + length(colors1);
    scatter(centers1(:,1),centers1(:,2),[],colors1)
    % viscircles(centers,repmat(5,size(centers,1),1),'LineWidth',0.5,'Color',colors1);
%     avgcenters1 = avgposc(avgposc(:,3)==k,1:2);
%     avgcolors1 = avgposc(avgposc(:,3)==k,4:6);
%     avgcenters(p:p+length(avgcenters1)-1,:) = avgcenters1;
%     avgcolors(p:p+length(avgcolors1)-1,:) = avgcolors1;
%     scatter(avgcenters1(:,1),avgcenters1(:,2),[],avgcolors1);
%     % viscircles(avgcenters,repmat(5,size(avgcenters,1),1),'LineWidth',0.5,'Color',avgcolors1);
%     p = p + length(avgcolors1);
    pause(0.5);
end