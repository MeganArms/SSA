function [C,B,F,M,R] = collate(objs_link)

% Collate trajectories
numTrajs = max(objs_link(6,:));
coords = objs_link(1:2,:);
brightness = objs_link(3,:);
% molInds = objs_link(4,:);
frames = objs_link(5,:);
trkID = objs_link(6,:);
indices = 1:length(objs_link);
C = cell(numTrajs,1); B = C; F = C; M = C;
R = struct;

for I = 1:numTrajs
    C{I} = coords(:,trkID == I); 
    B{I} = brightness(trkID == I);
    F{I} = frames(trkID == I);
    M{I} = indices(trkID == I);
    R(I).trajectory = indices(trkID == I);
end