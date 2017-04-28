function S = particleSize(obj,varargin)

% Analyze the input Molecule file for the frequency of fluorescence intensities 
% Input data should be the Molecule file.
% Output data is the histogram of intensities for this Molecule file, and
% the volume integrals of all the molecules found in this video. The
% volIntegrals variable should be concatenated over all the trials (of one
% experimental condition) to get the overall histogram of intensities.

Molecule = obj.Molecule;
if length(varagin) >= 1 && ~strcmp(varargin{1},'longTraj') || isempty(varargin)
    longTraj = obj.Result;
end
pixelSize = obj.Option.pixelSize;
R = obj.Option.spotR;
% Analyze molecules that appear on multiple frames and get their average
% intensities
MoleculeIndices = (1:length(Molecule))';
MoleculesAnalyzed = zeros(length(Molecule),1); 
S = zeros(2*length(longTraj),9);
L = length(longTraj);
for i = 1:length(longTraj)
    Size = zeros(length(longTraj(i).trajectory),2);
    Width = zeros(length(longTraj(i).trajectory),1);
    for j = 1:length(longTraj(i).trajectory)
        if j == 1
            mInitial = longTraj(i).trajectory(1);
        end
        if ~isnan(longTraj(i).trajectory(j))
           mIndex = longTraj(i).trajectory(j);
           if isfield(Molecule,'fit') %gaussian
               Width(j) = Molecule(mIndex).fit.sigma;
               if Width(j) == 0.1*R*pixelSize || Width(j) == 0.6*R*pixelSize % at the boundaries of the fit
                   continue
               end
               [Size(j,1), Size(j,2)] = fitVolume(mIndex,Molecule);
               obj.Molecule(mIndex).volume = Size(j,1);
               obj.Molecule(mIndex).maxInt = Size(j,2); % Change this to the sum of 9 central
           elseif isfield(Molecule,'area') %centroid -> fast
               a = Molecule(mIndex).area;
               Width(j) = sqrt(a/pi);
               Size(j,1) = Molecule(mIndex).volume;
               Size(j,2) = Molecule(mIndex).maxInt;
           end
           MoleculesAnalyzed = MoleculesAnalyzed + MoleculeIndices.*(MoleculeIndices == longTraj(i).trajectory(j));
        end
    end
    d = diag(squareform(pdist(Molecule(mInitial).Coords)),1);
    if ge(length(longTraj(i).trajectory),4)
        % S(i,:) = [max(Size(:,1)), max(Size(:,2)), length(longTraj(i).trajectory), mean(Width)];
        S(i,:) = [max(Size(:,1)), std(Size(:,1)), max(Width), std(Width), max(Size(:,2)), std(Size(:,2)), mean(d,'omitnan'), std(d,'omitnan'), length(longTraj(i).trajectory)];
    elseif length(longTraj(i).trajectory)==3
        S(i+L,:) = [max(Size(:,1)), std(Size(:,1)), max(Width), std(Width), max(Size(:,2)), std(Size(:,2)), mean(d,'omitnan'), std(d,'omitnan'), length(longTraj(i).trajectory)];
    end

end

%Remove all blank lines in the matrix (corresponding to molecules that
%appear in less than 3 frames)
for i=size(S,1):-1:1
    if all(S(i,:))==0
        S(i,:)=[];
    end
end

% Analyze molecules that appear on only one frame
% 
% MoleculesRemaining = MoleculeIndices(MoleculeIndices ~= MoleculesAnalyzed);
% starti = length(S)+1;
% S = [S; zeros(length(MoleculesRemaining),4)];
% for i = starti:length(S)
%     mIndex = MoleculesRemaining(i-starti+1);
%     if isfield(Molecule,'fit')
%         [volInt, maxInt] = fitVolume(mIndex, Molecule);
%         S(i,:) = [volInt, maxInt, 1, Molecule(mIndex).fit.sigma];
%         obj.Molecule(mIndex).volume = volInt;
%         obj.Molecule(mIndex).maxInt = maxInt;
%         obj.Molecule(mIndex).width = obj.Molecule(mIndex).fit.sigma;
%     elseif isfield(Molecule,'area')
%         volInt = Molecule(mIndex).volume;
%         maxInt = Molecule(mIndex).maxInt;
%         S(i,:) = [volInt, maxInt, 1, sqrt(Molecule(mIndex).area/pi)];
%     end
% end

% Analyze molecules that appear on at least 3 frames

% Scale the intensities to counts from counts*um^2
[obj.Intensity, ~] = VItransform(S, 190, 13);

    function [scaled, newVI] = VItransform(volIntegral, musigma, sdsigma)
        
        % Scale intensity integrals: Intensities were scaled over discrete bins of
        % 160 x 160 nm, so must be scaled down to a point, i.e. divided by pixelSize^2
        % nm^2. MEANSIGMA and SDSIGMA were manually determined by sampling
        % several image movie files, and this must done for each new
        % testing condition.
        
        scaled = volIntegral;
        scaled(:,1) = volIntegral(:,1)./pixelSize^2;
        scaled(:,2) = volIntegral(:,2)./pixelSize^2;
        
        % Create a new matrix with only the molecules that are +/- three standard deviation from
        % the mean sigma fit value.
        % varname = inputname(1);
        % musigma = 150; sdsigma = 90;
        minsigma = -3*sdsigma + musigma;
        maxsigma = 3*sdsigma + musigma;
        rows2keep = scaled(:,4) > minsigma & scaled(:,4) < maxsigma;
        newVI = scaled(rows2keep,:);
        
    end

    function [volumeInt, maxInt] = fitVolume(M, Molecule)
        
        % Find the volume of the Gaussian fit for molecule M using the fit
        % f = A*exp(-((x-x_0)^2 + (y-y_0)^2)/2/sigma^2) + z_0
        % The units of the parameters are in distance (usually microns) and
        % intensity levels
        
        A = Molecule(M).fit.A;
        sigma = Molecule(M).fit.sigma;
        
        % UPDATE so as to not include the background, z_0 level
        % syms x y A x_0 y_0 sigma z_0
        % f = A*exp(-((x-x_0)^2 + (y-y_0)^2)/2/sigma^2);
        % int_fx = int(f, x, [x_0-3*sigma, x_0 + 3*sigma]);
        % int_f = int(int_fx, y, [y_0 - 3*sigma, y_0 + 3*sigma])
        %
        % int_f =
        %
        % 2*A*pi*sigma^2*erf((3*2^(1/2))/2)^2
        
        volumeInt = 2*A*pi*sigma^2*erf((3*2^(1/2))/2)^2;
        maxInt = A;
    end

% Plot intensity versus sigma and hist of intensities if plotting is on
if length(varargin) >= 1 && strcmp(varargin{2},'on') || strcmp(varargin{1},'on')
    figure, 
    subplot(1,2,1), plot(obj.Intensity(:,4),obj.Intensity(:,1),'o');
    xlabel('Standard Deviation of Fits ({\mu}m)'), ylabel('Intensity (a.u.)');
    subplot(1,2,2), hist(obj.Intensity(:,1),100,'FaceColor','c');
    xlabel('Intensity (a.u.)'),ylabel('f(Intensity)');
    title(obj.filename(end-52:end-32));
end

% Get histogram of each fluorescence for each bin of visible time. Plot if
% plotting is on
% intensityDecayHist = hist3(volIntegrals,[100, 100]);
% if length(varargin) >= 1 && strcmp(varargin{1},'on')
%     v2plot = [volIntegrals(:,1)./length(volIntegrals), volIntegrals(:,2)*100];
%     figure, hist3(v2plot,[100, 100]);
%     xlabel('Fluorescence Intensity'), ylabel('Time Molecule Visible (ms)'), zlabel('Probability');
%     set(gcf,'renderer','opengl');
%     set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
% end

% % Get frequency (histogram) of each fluorescence intensity and plot, if 
% % plotting is on
% intensityHist = histc(volIntegrals(:,1), linspace(min(volIntegrals(:,1)),max(volIntegrals(:,1))));
% if length(varargin) >= 1 && strcmp(varargin{1},'on')
%     figure, bar(linspace(min(volIntegrals(:,1)),max(volIntegrals(:,1))), intensityHist, 'histc');
%     xlabel('Intensity Counts (a.u.)'); ylabel('Frequency of Intensity Counts');
% end
end
