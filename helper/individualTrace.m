function traces = individualTrace(obj,longTraj)

% Get the intensity traces of the connected trajectories. LongTraj should
% be the output of connectEvents. The difference between connectEvents and
% connectShortTraj is that there is no time limit on the connections
% between events in connectEvents, and the distance overwhich they might be
% connected is much larger than in connectShortTraj.

% Read video into 3D array for quick access
iminf = imfinfo(obj.filename);
N = length(iminf);
img = zeros(iminf(1).Width,iminf(1).Height,N);

for i = 1:N
    img(:,:,i) = imread(obj.filename,i);
end

% Obtain traces
n = length(longTraj);
traces = zeros(n,N);
for i = 1:n
    m = length(longTraj(i).trajectory);
    for j = 1:N
        if  j <= m && ~isnan(longTraj(i).trajectory(j))
            molInd = longTraj(i).trajectory(j);
            coords = obj.Molecule(molInd).coordinate;
            traces(i,j) = sum(sum(img(coords(1)-1:coords(1)+1,coords(2)-1:coords(2)+1,j)));
        else
            traces(i,j) = sum(sum(img(coords(1)-1:coords(1)+1,coords(2)-1:coords(2)+1,j)));
        end
        
    end
end
% traces = convert2double(traces);

