function [C,B,F,M,S,R,Frame] = collate(objs_link)

% Collate trajectories
numTrajs = max(objs_link(6,:));
coords = objs_link(1:2,:);
brightness = objs_link(3,:);
% molInds = objs_link(4,:);
frames = objs_link(5,:);
trkID = objs_link(6,:);
sigma = objs_link(7,:);
indices = 1:length(objs_link);
C = cell(numTrajs,1); B = C; F = C; M = C; S = C;
R = struct;
% f = @(x,y)x(y);
h = waitbar(0,'Collating data...');
for I = 1:numTrajs
    locs = trkID == I;
    C{I} = coords(:,locs); 
    B{I} = brightness(locs);
    F{I} = frames(locs);
    M{I} = indices(locs);
    S{I} = sigma(locs);
    R(I).trajectory = indices(locs);
    waitbar(I/numTrajs)
end
close(h)
clear h
Nframes = max(objs_link(5,:));
Frame = cell(Nframes,1); % The objs_link index of the molecules on each frame
for i = 1:Nframes
    Frame{i} = find(objs_link(5,:)==i);
end