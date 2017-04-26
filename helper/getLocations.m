function projectedCoordinates = getLocations(obj)

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
%       - C - 
Molecule = obj.Molecule;
Option = obj.Option;

projectedCoordinates = zeros(length(Molecule),1);
cs = Option.pixelSize/1000; % Coordinate scaling factor from pixels to microns
ps = 1/1000; % Parameter scaling factor from nanometers to microns

for i = 1:length(Molecule)
    if isfield(Molecule,'fit')
        projectedCoordinates(i,1) = cs*Molecule(i).coordinate(1) ...
            + ps*Molecule(i).fit.y0;
        projectedCoordinates(i,2) = cs*Molecule(i).coordinate(2) ...
            + ps*Molecule(i).fit.x0;
    elseif isfield(Molecule,'centroid')
        projectedCoordinates(i,1) = cs*Molecule(i).coordinate(1) ...
            + ps*Molecule(i).centroid(1);
        projectedCoordinates(i,2) = cs*Molecule(currentMolecule).coordinate(2) ...
            + ps*Molecule(currentMolecule).centroid(2);
    end
end
end