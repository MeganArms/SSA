function [spotEvents, eventfreq] = eventlinks01(tol,ttol,et,C,M,B,varargin)

% Organize inputs
if ~isempty(varargin)
	F = varargin{1};
	C = cellfun(@(x)x',C,'UniformOutput',false);
else
	F = cellfun(@(x)x(:,3),C,'UniformOutput',false);
end

%% Create searchable array
longCoords = C; 
clengths = zeros(length(C),1);
objnum = M;

% Columns of tosearch: Molecule number, last frame, i_last, j_last, first
% frame, i_first, j_first
tosearch = nan(length(longCoords),9);
for k = 1:length(longCoords)
    clengths(k) = size(longCoords{k},1);
    if clengths(k) <= 1
        longCoords{k} = [];
        objnum{k} = [];
    else
        tosearch(k,1) = k;
        tosearch(k,2) = objnum{k}(end);
        tosearch(k,3) = F{k}(end);
        tosearch(k,4) = C{k}(clengths(k),1);
        tosearch(k,5) = C{k}(clengths(k),2);
        tosearch(k,6) = objnum{k}(1);
        tosearch(k,7) = F{k}(1);
        tosearch(k,8) = C{k}(1,1);
        tosearch(k,9) = C{k}(1,2);
    end
end

%% Find connectable trajectories

ftol = ttol/et;
spotEvents = struct;
eventfreq = zeros(length(longCoords),1);
lastframes = unique(tosearch(:,3));
lastframes = lastframes(~isnan(lastframes));
toconnect = cell(length(longCoords),1);
for k = 1:length(longCoords)
    toconnect{k} = k;
end

h = waitbar(0,'Finding trajectories to connect');
for k = 1:length(lastframes)
    currlastframe = lastframes(k);
    rows_with_lastframe = tosearch(:,3) == currlastframe;
    seeds = tosearch(rows_with_lastframe,4:5);
    longCoords_inds = tosearch(rows_with_lastframe,1);
    in_tol_inds = find(tosearch(:,7) - currlastframe < ftol & tosearch(:,7) - currlastframe > 1 & ~rows_with_lastframe);

    [idx,D] = rangesearch(tosearch(in_tol_inds,8:9),seeds,tol);
    
    closest_traj = zeros(length(idx),1);
    closest_dist = zeros(length(idx),1);
    l = 1;
    while l <= length(idx)
        if ~isempty(idx{l})
            closest_traj(l) = idx{l}(1);
            closest_dist(l) = D{l}(1);
            l = l + 1;
        else
            longCoords_inds(l) = [];
            idx(l) = [];
            D(l) = [];
        end
    end
    closest_dist = closest_dist(closest_dist ~= 0);
    closest_traj = closest_traj(closest_traj ~= 0);
    
    for l = 1:length(idx)
        if ~isempty(idx{l})
%             toconnect{longCoords_inds(l)} = [toconnect{longCoords_inds(l)},in_tol_inds(closest_traj(l))];
            %         reveal_repeats = closest_traj(:) == closest_traj(l);
            s = 1;
            while s <= length(idx{l})
                closest = idx{l}(s);
                find_matches = closest_traj(:) == closest;
                if sum(find_matches) > 1
                    [~,min_ind] = min(closest_dist(find_matches));
                    A = find(find_matches);
                    idx_ind_of_traj = A(min_ind); % Index of trajectory to be extended in IDX
                    ind_of_traj = in_tol_inds(idx{idx_ind_of_traj}(1)); % Index of extending trajectory in LONGCOORDS
                    
                    % Remove other mathcing trajectories
                    find_matches(idx_ind_of_traj) = false;
                    A = find(find_matches);
                    for p = 1:length(A)
                        if length(idx{A}) > 1
                            where_other_close_trajs = idx{A(p)} ~= closest_traj(A(p));
                            other_close_trajs = idx{A(p)}(where_other_close_trajs);
                            closest_traj(A(p)) = other_close_trajs(1);
                            other_close_dist = D{A(p)}(where_other_close_trajs);
                            closest_dist(A(p)) = other_close_dist(1);
                        else
                            closest_dist(A) = NaN;
                            closest_traj(A) = NaN;
                            longCoords_inds(A) = NaN;
                            idx{A} = [];
                        end
                    end
                    
                    A = longCoords_inds(idx_ind_of_traj);
                    toconnect{A} = [toconnect{A},ind_of_traj];
                    tosearch(A,3:5) = tosearch(ind_of_traj,3:5);
                    tosearch(ind_of_traj,:) = nan(1,9);
                    s = s + 1;
                else
                    A = longCoords_inds(l);
                    I = in_tol_inds(closest);
                    toconnect{A} = [toconnect{A},I];
                    tosearch(A,3:5) = tosearch(I,3:5);
                    tosearch(I,:) = nan(1,9);
                    break
                end
            end
        end
    end
    waitbar(k/length(lastframes),h);
end
close(h)
clear h
%% Connect connectable trajectories

p = 1;
for k = 1:length(toconnect)
    if ~isempty(longCoords{k})
        traj2connect = toconnect{k};
        if length(traj2connect) > 1
            tmptraj = M{traj2connect(1)};
            for l = 1:length(traj2connect)-1
                currlastframe = F{traj2connect(l)}(end);
                nextfirstframe = F{traj2connect(l+1)}(1);
                spacer = nan(1,nextfirstframe - currlastframe - 1);
                tmptraj = [tmptraj,spacer,M{traj2connect(l+1)}];
            end
            spotEvents(p).trajectory = tmptraj;
        else
            spotEvents(p).trajectory = M{traj2connect};
        end
        spotEvents(p).coordinates = cat(1,longCoords{traj2connect});
        spotEvents(p).brightness = cat(2,B{traj2connect});
        spotEvents(p).frames = cat(2,F{traj2connect});
        spotEvents(p).eventfreq = length(traj2connect);
        eventfreq(p) = length(traj2connect);
        p = p + 1;
    end
end
eventfreq = eventfreq(eventfreq~=0);

end