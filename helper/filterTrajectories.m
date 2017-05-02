function [Result, censor] = filterTrajectories(obj,spots,exclude,th,d0,dstart)

% This function converts the molecular indices of the trajectories in
% obj.Result to trajectories with coordinates in microns, with [i,j] =
% [1,1] as the origin.
% 
%   Input: 
%       - OBJ - FSMIA object
%       - SPOTS - connected trajectories vector, output of connectEvents.m
%       to connect proteins past blinking.
%       - EXCLUDE - 'yes' will exclude trajectories with molecules
%       appearing on the first and/or last frame from any further analysis
%       with "coords" field of OBJ struct. 'no' will include all
%       trajectories.
%       - TH - Minimum number of frames that the object must appear on the surface
%       - D0 - Experiment start time in the format 'yyyyMMdd HHmmss'
%       - DSTART - Observation start time in the format 'yyyyMMdd HHmmss'
%   Output:
%       - RESULT - filtered trajectories - censored or padded as stipulated in
%       the input
%       - CENSOR - vector with censoring of trajectories. Censoring is "1" if
%       the molecule appears on the first frame, 0 otherwise.

t0 = datetime(d0,'InputFormat','yyyyMMdd HHmmss');
tstart = datetime(dstart,'InputFormat','yyyyMMdd HHmmss');
D = datenum(tstart - t0)*24*60*60;
incubation = floor(D/et);

Molecule = obj.Molecule;
Option = obj.Option;
L = length(spots);
T = length(obj.Frame);
if isfield(Molecule,'connectEvents')
    Result = struct;
    if strcmp(exclude,'yes') % When there is photobleaching detection, this case will also need a cens vector
        k = 1;
        for t = 1:length(spots)
            mStart = spots(t).trajectory(1);
            mEnd = spots(t).trajectory(end);
            if Molecule(mStart).frame ~= 1 && Molecule(mEnd).frame ~= T && length(spots(t).trajectory) >= th
                Result(k).trajectory = spots(t).trajectory;
                k = k + 1;
            end
        end
    elseif strcmp(exclude,'no')
        k = 1; cens = zeros(L,1);
        for t = 1:length(spots)
            mStart = spots(t).trajectory(1);
            mEnd = spots(t).trajectory(end);
            if length(spots(t).trajectory) >= th
                if Molecule(mEnd).frame == T
                    cens(k) = 1;
                elseif Molecule(mEnd).frame ~= T
                    cens(k) = 0;
                end
                Result(k).trajectory = spots(t).trajectory;
                if Molecule(mStart).frame == 1
                    Result(k).experimental = padarray(spots(t).trajectory,[0 incubation],'pre');
                else
                    Result(k).experimental = spots(t).trajectory;
                end
                k = k + 1;
            end
        end
        if nargout > 1
            censor = logical(cens);
        end
    end
else
    Result = struct;
    if strcmp(exclude,'yes')
        k = 1; t = 1;
        while k < length(Result)
            mindex = obj.Result(k).trajectory(1);
            if Molecule(mindex).frame ~= 1 && length(obj.Result(k).trajectory) >= th
                Result(t).trajectory = obj.Result(k).trajectory;
                t = t + 1;
            end
            k = k + 1;
        end
    elseif strcmp(exclude,'no')
        Result = obj.Result;
        cens = zeros(length(Result),1);
        for k = 1:length(Result)
            mindex = Result(k).trajectory(1);
            if Molecule(mindex).frame == T
                cens(k) = 1;
            elseif Molecule(mindex).frame ~= T
                cens(k) = 0;
            end
            if Molecule(mStart).frame == 1
                Result(k).experimental = padarray(spots(t).trajectory,[0 incubation],'pre');
            else
                Result(k).experimental = spots(t).trajectory;
            end
        end
        if nargout > 1
            censor = logical(cens);
        end
    end
end 

end