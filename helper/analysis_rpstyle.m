function driftobj = analysis_rpstyle(img)

threshopt = 1; thresh = 0.999; nsize = 7; try1pernhood = false; dimg = [];
[y, x] = calcthreshpts(img, threshopt, thresh, nsize, try1pernhood, dimg);
orientationstr = 'momentcalc';

% Get rid of maxima too close to the edge
lenx = 2*floor(nsize/2) + 1;  % 'floor' isn't really necessary, but this
% is the size of "nhood = getnhood(ste);" for a disk structuring
% element of that size.  Note that lenx is forced to be odd.
leny = lenx;  % in principle, could make different
edgeind = ((x < lenx/2) | (x > (size(img,2) - lenx/2))) | ...
    ((y < leny/2) | (y > (size(img,1) - leny/2)));
x(edgeind) = [];
y(edgeind) = [];

%% Determine neighborhoods
rect = zeros(length(x), 4);
% Compute the first local neighborhood to know the image size, and then all
% the rest.  A bit inelegant; could calculate all rect's at once...
% Skip if there are no objects to find
if ~isempty(x)
    rect(1,:) = [(round(x(1)) - floor(lenx/2)) (round(y(1)) - floor(leny/2)) (lenx-1) (leny-1)];
    % cropimg1 = imcrop(img, rect(1,:));
    cropimg1 = img(rect(1,2):(rect(1,2)+lenx-1), rect(1,1):(rect(1,1)+lenx-1));
    % all the other neighborhoods
    cropimg = repmat(cropimg1, [1 1 length(x)]); % to allocate memory
    for k = 2:length(x)
        rect(k,:) = [(round(x(k)) - floor(lenx/2)) (round(y(k)) - floor(leny/2)) (lenx-1) (leny-1)];
        % cropimg(:,:,k) = imcrop(img, rect(k,:));
        cropimg(:,:,k) = img(rect(k,2):(rect(k,2)+lenx-1), rect(k,1):(rect(k,1)+lenx-1));
    end
end
%% Compute "masses" (particle brightness)
% Calculate "mass" (intensity) in each neighborhood
% subtract background, estimated from mean of non-neighborhood px.
nhood = getnhood(strel('disk', floor(nsize/2),0));  % somewhat silly; 
        % could avoid call to image processing toolbox if needed
        % Note e.g. nhood = getnhood(strel('disk', floor(5/2),0)) is
%              0     0     1     0     0
%              0     1     1     1     0
%              1     1     1     1     1
%              0     1     1     1     0
%              0     0     1     0     0
meanbkg = zeros(1, length(x));
savemass = zeros(1, length(x));
for k = 1:length(x)
    tempreg = cropimg(:,:,k);
    cropreg = tempreg(nhood);
    bkgreg  = tempreg(~nhood);
    meanbkg(k) = mean(bkgreg(:));
    savemass(k) = sum(cropreg(:)) - sum(nhood(:))*meanbkg(k);
end
%% Refine positions
% Do refinement (find center) around each local max
% If there are many points, use radialcenter_stk.m rather than radialcenter.m
% for radial symmetry method -- even faster!
xcent = zeros(1,length(x));
ycent = zeros(1,length(x));
sigma = zeros(1,length(x));
lsumx = 1:size(cropimg,2);
lsumy = 1:size(cropimg,1);
Lx = lsumx(end);
Ly = lsumy(end);
meand2 = zeros(1, length(x));  % this will only be filled in for 'radial' localization
% Radial-symmetry based fit -- fast, accurate
% If <10 points, use radialcenter_stk.m for extra speed (avoid
% redundant grid calculations); else radialcenter.m
if length(x) < 10
    for j = 1:length(x)
        [xcent(j), ycent(j), sigma(j), meand2(j)] = radialcenter(cropimg(:,:,j));
    end
else
    [xcent, ycent, sigma, meand2] = radialcenter_stk(cropimg) ;
end
% Is the center within reasonable bounds?
% If not, replace with centroid
% frequency of bad cases ~ 1/100,000 !  (Extremely rare)
% See notes Oct. 26, 2011
% This conditional statement can slow things; Delete?
badcase = find(abs(xcent - Lx/2)>1.5*Lx | abs(ycent - Ly/2)>1.5*Ly);
for j = badcase
    ci = cropimg(:,:,j);
    xcent(j) = sum(sum(ci) .* lsumx) / sum(sum(ci));
    ycent(j) = sum(sum(ci,2) .* lsumy') / sum(sum(ci));
end
    
% center position relative to image boundary
% Note that the *center* of the upper left pixel is (1,1)
xn = xcent + rect(:,1)' - 1; % -1 is to correct for matlab indexing
yn = ycent + rect(:,2)' - 1;

%% If desired, calculate the ellipsoid orientation
if ~isempty(orientationstr)
    theta = zeros(1,length(x));
    ra = zeros(1,length(x));
    rb = zeros(1,length(x));
    switch lower(orientationstr)
        % Note that 'none' has already been turned into an empty
        % orientation string
        case {'momentcalc'}
            % simple moment calculation
            % Redundant grid calculations, etc., that could be sped up...
            for j = 1:length(x)
                [~, theta(j), r] = simpleellipsefit(cropimg(:,:,j), [xcent(j) ycent(j)], false, true);
                ra(j) = r(1);
                rb(j) = r(2);
            end
        otherwise
            errordlg('Unknown orientation method! [fo5_rp.m]');
    end
end

%% Create objs matrix
nrows = 8;
driftobj = zeros(nrows, length(xn));
driftobj(1,:) = xn;
driftobj(2,:) = yn;
driftobj(3,:) = savemass;
driftobj(4,:) = 1:length(x);
driftobj(7,:) = sigma;
driftobj(8,:) = meand2;  % zero for all methods other than 'radial'
driftobj(9,:) = theta;
driftobj(10,:) = ra;
driftobj(11,:) = rb;