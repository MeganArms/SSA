function [Result, cens] = filterTrajectories_rp(objs_link,R,exclude,th,d0,dstart,et,Nframes)

% This function converts the molecular indices of the trajectories in
% obj.Result to trajectories with coordinates in microns, with [i,j] =
% [1,1] as the origin.
% 
%   Input: 
%       - OBJ_LINKS - output of TrackingGUI_rp
%       - R - connected trajectories vector, output of connectEvents.m
%       to connect proteins past blinking.
%       - EXCLUDE - 'yes' will exclude trajectories with molecules
%       appearing on the first and/or last frame from any further analysis
%       with "coords" field of OBJ struct. 'no' will include all
%       trajectories. If this is used with robust survival analysis, set to 'no'.
%       - TH - Minimum number of frames that the object must appear on the surface
%       - D0 - Experiment start time in the format 'yyyyMMdd HHmmss'
%       - DSTART - Observation start time in the format 'yyyyMMdd HHmmss'
%       - ET - Time between the start of exposures in seconds
%   Output:
%       - RESULT - filtered trajectories - censored or padded as stipulated in
%       the input
%       - CENSOR - vector with censoring of trajectories. Censoring is "1" if
%       the molecule appears on the first frame, 0 otherwise.

t0 = datetime(d0,'InputFormat','yyyyMMdd HHmmss');
tstart = datetime(dstart,'InputFormat','yyyyMMdd HHmmss');
D = datenum(tstart - t0)*24*60*60;
incubation = floor(D/et);

T = Nframes;
Result = struct;

if strcmp(exclude,'yes')
    k = 1; cens = zeros(length(R),1);
    for t = 1:length(R)
        mindex = R(t).trajectory(1);
        if objs_link(5,mindex) ~= 1 && length(R(t).trajectory) >= th
            Result(k).trajectory = R(t).trajectory;
            if objs_link(5,mindex) == T
                cens(k) = 1;
            elseif objs_link(5,mindex) ~= T
                cens(k) = 0;
            end
            k = k + 1;
        end
    end
elseif strcmp(exclude,'no')
    Result = struct;
    cens = zeros(length(R),1); k = 1;
    for t = 1:length(R)
        mStart = R(t).trajectory(1);
        mEnd = R(t).trajectory(end);
        if length(R(t).trajectory) >= th
            Result(k).trajectory = R(t).trajectory;
            if objs_link(5,mEnd) == T
                cens(k) = 1;
            elseif objs_link(5,mEnd) ~= T
                cens(k) = 0;
            end
            if objs_link(5,mStart) == 1
                Result(k).experimental = padarray(R(t).trajectory,[0 incubation],'pre');
            else
                Result(k).experimental = R(t).trajectory;
            end
            k = k + 1;
        end

    end
end
if nargout > 1
    cens = cens(1:length(Result));
end


end