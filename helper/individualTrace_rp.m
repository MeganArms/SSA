function traces = individualTrace_rp(filename,longTraj,th,objs_link)

% Get the intensity traces of the connected trajectories. LongTraj should
% be the output of connectEvents (spotEvents). The difference between connectEvents and
% connectShortTraj is that there is no time limit on the connections
% between events in connectEvents, and the distance overwhich they might be
% connected is much larger than in connectShortTraj.
% th is the minimum trajectory length

% Read video into 3D array for quick access
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

% [m,n,N] = size(img);
% Obtain traces
L = length(longTraj);
traces = zeros(L,N);
for k = 1:L
    TOI = longTraj(k).trajectory;
    t = length(TOI);
    startFr = objs_link(5,find(~isnan(TOI),1));
    if t >= th
        if startFr == 1
            for l = 1:N
                if l <= t && ~isnan(longTraj(k).trajectory(l))
                    coords = floor(objs_link(1:2,longTraj(k).trajectory(l)));
                    traces(k,l) = objs_link(3,longTraj(k).trajectory(l))+1000;
                else
                    ind = sub2ind([m n N],coords(1),coords(2),l);
                    traces(k,l) = img(ind);
                end
            end
        else
            coords = floor(objs_link(1:2,longTraj(k).trajectory(1)));
            for l = 1:startFr-1
                ind = sub2ind([m n N],coords(1),coords(2),l);
                traces(k,l) = img(ind);
            end
        end
    end
end
traces = traces(sum(traces,2)>0,:);
% traces = convert2double(traces);

