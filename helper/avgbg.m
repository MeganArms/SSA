function I = avgbg(varargin)

% AVGBG Gets the average background intensity of each frame of a video
% file. Any objects in the foreground and removed and the average intensity
% of background is taken. 
% INPUT:
%   VARARGIN - 
%       - If there is no input, user will be prompted to select the image
%       file
%       - FILENAME is accepted here as full path to image file
% OUTPUT:
%   I - vector of average background intensities on each frame

if isempty(varargin)
    [filename,user_canceled] = imgetfile;
    if user_canceled == 1
        I = NaN;
        disp('No image selected')
        return
    end
else
    filename = varargin{1};
end

Nframes = length(imfinfo(filename));
avgint = cell(Nframes,2);
h = waitbar(0,'Analyzing background levels...');
for i = 1:Nframes
    img = imread(filename,i);
    I1 = imopen(img,strel('diamond',11));
    avgint{i,1} = histcounts(I1(:),'Normalization','pdf');
    avgint{i,2} = mean(I1(:));
    waitbar(i/Nframes)
end
close(h)
I = cell2mat(avgint(:,2));