function traces = individualTrace_rp(filename,longTraj,C,th)

% Get the intensity traces of the connected trajectories. LongTraj should
% be the output of connectEvents (spotEvents). The difference between connectEvents and
% connectShortTraj is that there is no time limit on the connections
% between events in connectEvents, and the distance overwhich they might be
% connected is much larger than in connectShortTraj.
% th is the minimum trajectory length

% Read video into 3D array for quick access
iminf = imfinfo(filename);
N = length(iminf);
img = zeros(iminf(1).Width,iminf(1).Height,N);
C = C(cellfun(@length,C) >= th);

for i = 1:N
    img(:,:,i) = imread(filename,i);
end

% Obtain traces
n = length(longTraj);
traces = zeros(n,N);
for i = 1:n
    m = length(longTraj(i).trajectory);
    for j = 1:N
        if  j <= m && ~isnan(longTraj(i).trajectory(j))
            coords = floor(C{i}(:,j));
            traces(i,j) = img(coords(1),coords(2),j);
        else
            traces(i,j) = img(coords(1),coords(2),j);
        end
        
    end
end
% traces = convert2double(traces);

