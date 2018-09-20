function [survFuncs,intensities,numbers] = generateoutput(et,Nframes,PBflag,ttol,varargin) 
% GENERATEOUTPUT
% Inputs - 
%   ET: The time between the start of individual frames. Either the
%   exposure time or interval time.
%   NFRAMES: The number of frames in each video.
%   PBFLAG: Set to FALSE for experimental adsorption videos. Set to TRUE
%   for videos of chemisorption or other photobleaching control videos.
%   This setting only uses the particles that appear on the first frame for
%   the analysis, assuming that these will be fully adsorbed.
%   TTOL: The maximum time allowed between bright events for connection to
%   occur across dark events.
%   VARARGIN:
%       Variable argument inputs, requiring all of the below variables, but
%       can be repeated for each video trial of the same condition.
%       Example: 
%       [survFuncs,intensities,numbers] = generateoutput(et,Nframes,false,5...
%           C1,F1,objs_link1,spotEvents1,eventfreq1,...
%           C2,F2,objs_link2,spotEvents2,eventfreq2);
%
%       C: Cell array of paired coordinates of particles separated by 
%       trajectory.
%       F: Cell array of frames of particles separated by trajectory.
%       OBJS_LINK: key output from TRACKING_GUIRP analysis that contains
%       parameters describing the detected objects.
%       SPOTEVENTS: Particle trajectories that have been connected across
%       dark events. Output of EVENTLINKS9 or SEGMENTTRAJ.
%       EVENTFREQ: The frequency of events per spot.
%       Output of EVENTLINKS9. Perhaps archaic calculation method.

% Sort inputs
if iscell(varargin{end})
    survFuncs = varargin{end-2};
    intensities = varargin{end-1};
    numbers = varargin{end};
    C = cat(1,varargin{1:5:end-3});
    F = cat(1,varargin{2:5:end-3});
    objs = varargin(3:5:end-3);
    spots = varargin(4:5:end-3);
    eventfreq = cat(1,varargin{5:5:end-3});
else
    C = cat(1,varargin{1:5:end});
    F = cat(1,varargin{2:5:end});
    objs = varargin(3:5:end);
    spots = varargin(4:5:end);
    eventfreq = cat(1,varargin{5:5:end});
end

% Generate output from helper functions
[sfc,tc,semc,bc,numec,~,~] = crtd(Nframes,et,PBflag,spots,objs);
[fc,xc,reside,cens,floc,fupc] = km(Nframes,et,PBflag,spots);
[dtime,firstframe,lastframe,brightness,avgb,rangeb,maxb,rcdfb,eventrate] = factors(Nframes,et,PBflag,spots,objs);
[f,x,~,~,~,CC,b,xx,nume] = residenceTime(C,F,Nframes,et);
[dtime,edges] = histcounts(dtime,logspace(-2,log(ttol),100),'Normalization','pdf');
x_d = edges(1:end-1);

num1 = sum(firstframe == 1); 
fffract = num1/length(firstframe);
numlast = sum(lastframe == Nframes);
lffract = numlast/length(lastframe);
numtrajc = sum(cellfun(@length,spots));
numobj = sum(cellfun(@(x)size(x,2),objs));
numtraj = length(C);
xc(1) = 0; x(1) = 0;

numperframe = histcounts([F{:}],1:Nframes+1);
frames = [];
for k = 1:length(spots)
    spotE = spots{k};
    frames = [frames,spotE(:).frames];
end
numperframec = histcounts(frames,1:Nframes+1);

% Collect output into cell arrays
if ~iscell(varargin{end})
    survFuncs = {sfc,tc,semc,bc,numec,fc,xc,reside,cens,CC,xx,b,nume,f,x,floc,fupc};
%     intensities = {icounts(:,1),icounts(:,2),icounts(:,3),icounts(:,4),x_d,dtime,avgb,rangeb,maxb,rcdfb,brightness};
    intensities = {x_d,dtime,avgb,rangeb,maxb,rcdfb,brightness};
    numbers = {firstframe,lastframe,num1,fffract,numlast,lffract,numtrajc,numtraj,numobj,eventfreq,Nframes,numperframe,numperframec,eventrate};
else
    newIndex = size(survFuncs,1)+1;
    survFuncs(newIndex,:) = {sfc,tc,semc,bc,numec,fc,xc,reside,cens,CC,xx,b,nume,f,x,floc,fupc};
%     intensities(newIndex,:) = {icounts(:,1),icounts(:,2),icounts(:,3),icounts(:,4),x_d,dtime,avgb,rangeb,maxb,rcdfb,brightness};
    intensities(newIndex,:) = {x_d,dtime,avgb,rangeb,maxb,rcdfb,brightness};
    numbers(newIndex,:) = {firstframe,lastframe,num1,fffract,numlast,lffract,numtrajc,numtraj,numobj,eventfreq,Nframes,numperframe,numperframec,eventrate};
end

end