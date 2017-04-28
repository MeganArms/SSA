function FineScan(obj,RawImage,k)
% FINESCAN(MOLPIXELIDX,RAWIMAGE) gets the detailed information for molecules
% identified in RAWIMAGE

Option = obj.Option;
R = Option.spotR;   % radius (pixel) of diffraction limited spot
[molPixelIdx,BW] = RoughScan(obj,RawImage,k);
img = double(RawImage);
obj.Option.bg = min(img(:));%-std(img(:));
NumMolecule = length(obj.Molecule);

for k = 1:length(molPixelIdx)
    if isempty(molPixelIdx{k})
        break
    end
    i = molPixelIdx{k}(1);
    j = molPixelIdx{k}(2);
    subImage = img(i-R:i+R,j-R:j+R);
    BW_sub = BW(i-R:i+R,j-R:j+R);
    CC_sub = bwconncomp(BW_sub); %output struct 
    
    % Deal with the above threshold pixels in the peripheral of subimage
    N = numel(CC_sub.PixelIdxList); % number of objects in a subimage - should be 1
    if N > 1
        center_idx = 2*R^2+2*R+1;
        for l = 1:N
            pixIdxList = CC_sub.PixelIdxList{l};
            if ~ismember(center_idx,pixIdxList)
                subImage(pixIdxList) = min(min(subImage));
            end
        end
    end
    if strcmp(Option.fitting,'fast')
        % Perform centroid fitting. Subtract the location of the center
        % pixel to convert it to the distance from the center of the pixel.
        % Eliminate potential molecules in ROI
        edgeThreshold = Option.threshold;
        edgeImage = subImage; edgeImage(3:end-2, 3:end-2) = 0;
        subImage(edgeImage > edgeThreshold) = obj.Option.bg;
        centroid = regionprops(true(size(subImage)),subImage,'WeightedCentroid');
        s = centroid.WeightedCentroid(2)-R-1+0.5;
        t = centroid.WeightedCentroid(1)-R-1+0.5;
        obj.Molecule(NumMolecule+k).centroid = [s,t]*obj.Option.pixelSize; %centroid in nm
        if N > 1
            lengths = zeros(1,N);
            for l = 1:N
                lengths(l) = length(CC_sub.PixelIdxList);
            end
            [~,maxidx] = max(lengths);
            pxlist = CC_sub.PixelIdxList{maxidx};
        else
            pxlist = CC_sub.PixelIdxList{1};
        end
        obj.Molecule(NumMolecule+k).volume = sum(sum(subImage(pxlist) - min(min(subImage))))*max(img(:));    
        obj.Molecule(NumMolecule+k).area = length(pxlist)*Option.pixelSize^2;
        obj.Molecule(NumMolecule+k).maxInt = max(max(subImage))*max(img(:));
    elseif strcmp(Option.fitting,'slow')
        % Eliminate potential molecules in ROI
        edgeThreshold = Option.threshold;
        edgeImage = subImage; edgeImage(3:end-2, 3:end-2) = 0;
        subImage(edgeImage > edgeThreshold) = obj.Option.bg;
         try
             [obj.Molecule(NumMolecule+k).fit,obj.Molecule(NumMolecule+k).gof] = fit2D(obj,subImage);
         catch %if error in "try"
            disp('Unable to fit 2D Gaussian for the following molecule with adjusted ROI:');
            fprintf('%d, %d\n',i,j);
            disp(subImage);
            continue
        end

    end
    obj.Molecule(NumMolecule+k).coordinate = [i j];
end

end