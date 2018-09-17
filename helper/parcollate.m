function [C,B,F,M,S,R,Frame] = parcollate(objs_link)

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
% h = waitbar(0,'Collating data...');
parfor I = 1:numTrajs
    locs = trkID == I;
    coordstmp = coords;
    brightnesstmp = brightness;
    framestmp = frames;
    indicestmp = indices;
    sigmatmp = sigma;
    
%     for k = 1:length(locs)
%         l = locs(k);
%         c(1,k) = coords(1,l);
%         c(2,k) = coords(2,l);
%         b(k) = brightness(l);
%         f(k) = frames(l);
%         m(k) = indices(l);
%         s(k) = sigma(l)
%         r(k) = indices(l);
%     end
% 
        c = coordstmp(:,locs);
        b = brightnesstmp(locs);
        f = framestmp(locs);
        m = indicestmp(locs);
        s = sigmatmp(locs);
        r = indicestmp(locs);   
   

    C{I} = c; 
    B{I} = b;
    F{I} = f;
    M{I} = m;
    S{I} = s;
    R(I).trajectory = r;
%     waitbar(I/numTrajs)
end
% close(h)
% clear h
Nframes = max(objs_link(5,:));
Frame = cell(Nframes,1); % The objs_link index of the molecules on each frame
for i = 1:Nframes
    Frame{i} = find(objs_link(5,:)==i);
end