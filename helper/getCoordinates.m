function [C, cens] = getCoordinates(obj,spots,exclude)

% This function converts the molecular indices of the trajectories in
% obj.Result to trajectories with coordinates in microns, with [i,j] =
% [1,1] as the origin.
% 
%   Input: 
%       - OBJ - FSMIA object
%       - EXCLUDE - 'yes' will exclude trajectories with molecules
%       appearing on the first and/or last frame from any further analysis
%       with "coords" field of OBJ struct. 'no' will include all
%       trajectories.
%   Output:
%       - C - cell array in microns. Additionally stored with the Molecule
%       whose index is first in the trajectory.Each cell contains a matrix
%       with the I and J coordinates of the points in the trajectory.
%       Column 1 contains I, and column 2 contains J. For example, the 17th
%       path is located in the 17th cell. The coordinates of this path are
%       accessed as gaussCoordinates{17} if you want the coordinate values
%       (e.g. [371, 124; 371, 123]), or as gaussCoordinates(17) if you want
%       the matrix containing the coordinate values (e.g. [2x2 double]).
%       - CENS - vector with censoring of trajectories. Censoring is "1" if
%       the molecule appears on the first frame, 0 otherwise.

Molecule = obj.Molecule;
Option = obj.Option;
L = length(spots);
if isfield(Molecule,'connectEvents')
    Result = struct;
    if strcmp(exclude,'yes')
        k = 1;
        for i = 1:length(Molecule)
            if ~isempty(Molecule(i).connectEvents) && Molecule(i).frame ~= 1 && Molecule(Molecule(i).connectEvents(end)).frame ~= length(obj.Frame)
                Result(k).trajectory = Molecule(i).connectEvents;
                k = k + 1;
            end
        end
    elseif strcmp(exclude,'no')
        k = 1; cens = zeros(L,1);
        for i = 1:length(Molecule)
            if ~isempty(Molecule(i).connectEvents)
                if Molecule(i).frame == 1
                    cens(k) = 1;
                elseif Molecule(i).frame ~= 1
                    cens(k) = 0;
                end
                Result(k).trajectory = Molecule(i).connectEvents;
                k = k + 1;
            end
        end
    end
else
    Result = obj.Result;
    cens = zeros(length(Result),1);
    for k = 1:length(Result)
        if Molecule(Result(k).trajectory(1)).frame == 1
            cens(k) = 1
        elseif Molecule(Result(k).trajectory(1)).frame ~= 1
            cens(k) = 0;
        end
    end
end 
cens = logical(cens);
C = cell(length(Result),1);
cs = Option.pixelSize/1000; % Coordinate scaling factor from pixels to microns
ps = 1/1000; % Parameter scaling factor from nanometers to microns

for i = 1:length(Result)
    moleculesInTrajectory = Result(i).trajectory; % Get molecular indices
    moleculeTrajectoryCoordinates = zeros(length(moleculesInTrajectory),2); % allocate array for the number of frames that the trajectory crosses
    for j = 1:length(moleculesInTrajectory)
        currentMolecule = moleculesInTrajectory(j); % Progress through each of the molecular indices
        if isnan(currentMolecule)
            currentMolecule = previousMolecule; % If the current molecule is too dim, use the coordinates of the previous molecule
        else
            previousMolecule = currentMolecule; % If the current molecule is bright enough, save it as the previous molecule in case it goes dim in the next frame
        end        
        % Save the precise coordinate for each of the molecules
        if isfield(Molecule,'fit')
            moleculeTrajectoryCoordinates(j,1) = cs*Molecule(currentMolecule).coordinate(1) ...
                + ps*Molecule(currentMolecule).fit.y0;
            moleculeTrajectoryCoordinates(j,2) = cs*Molecule(currentMolecule).coordinate(2) ...
                + ps*Molecule(currentMolecule).fit.x0;
        elseif isfield(Molecule,'centroid')
            moleculeTrajectoryCoordinates(j,1) = cs*Molecule(currentMolecule).coordinate(1) ...
                + ps*Molecule(currentMolecule).centroid(1);
            moleculeTrajectoryCoordinates(j,2) = cs*Molecule(currentMolecule).coordinate(2) ...
                + ps*Molecule(currentMolecule).centroid(2);
        end
    end
    C{i} = moleculeTrajectoryCoordinates;
    obj.Molecule(moleculesInTrajectory(1)).Coords = moleculeTrajectoryCoordinates;
end

end