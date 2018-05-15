function [survFuncs,intensities,numbers] = generateoutput(et,Nframes,C1,C2,F1,F2,spotEvents1,spotEvents2,objs_link1,objs_link2,eventfreq1,eventfreq2,varargin)

[sfc,tc,semc,bc,numec,~,~] = crtd2(spotEvents1,spotEvents2,Nframes,et,objs_link1,objs_link2);
[icounts,dtime,firstframe,lastframe,brightness,avgb,rangeb,maxb,rcdfb] = factors(spotEvents1,spotEvents2,Nframes,objs_link1,objs_link2);
[f,x,~,~,~,CC,b,xx,nume] = residenceTime([C1;C2],[F1;F2],Nframes,et);
[fc,xc,reside,cens] = km2(spotEvents1,spotEvents2,Nframes,et,objs_link1,objs_link2);
num1 = sum(firstframe == 1); 
fffract = num1/length(firstframe);
numlast = sum(lastframe == Nframes);
lffract = numlast/length(lastframe);
numtrajc = length(spotEvents1)+length(spotEvents2);
numobj = size(objs_link1,2)+size(objs_link2,2);
numtraj = length(C1)+length(C2);
xc(1) = 0; x(1) = 0;

if isempty(varargin)
    survFuncs = {sfc,tc,semc,bc,numec,fc,xc,reside,cens,CC,xx,b,nume,f,x};
    intensities = {icounts(:,1),icounts(:,2),icounts(:,3),icounts(:,4),dtime,avgb,rangeb,maxb,rcdfb,brightness};
    numbers = {firstframe,lastframe,num1,fffract,numlast,lffract,numtrajc,numtraj,numobj,[eventfreq1;eventfreq2],Nframes};
else
    survFuncs = varargin{1};
    newIndex = size(survFuncs,1)+1;
    intensities = varargin{2};
    numbers = varargin{3};
    
    survFuncs(newIndex,:) = {sfc,tc,semc,bc,numec,fc,xc,reside,cens,CC,xx,b,nume,f,x};
    intensities(newIndex,:) = {icounts(:,1),icounts(:,2),icounts(:,3),icounts(:,4),dtime,avgb,rangeb,maxb,rcdfb,brightness};
    numbers(newIndex,:) = {firstframe,lastframe,num1,fffract,numlast,lffract,numtrajc,numtraj,numobj,[eventfreq1;eventfreq2],Nframes};
end

% objs_link1 = objs_link;
% clear objs_link objs
% [C1,B1,F1,M1,S1,R1,Frame1] = collate(objs_link1);
% [spotEvents1, eventfreq1] = eventlinks9(2,C1,M1,B1,F1);
% 
% objs_link2 = objs_link;
% clear objs objs_link
% [C2,B2,F2,M2,S2,R2,Frame2] = collate(objs_link2);
% [spotEvents2, eventfreq2] = eventlinks9(2,C2,M2,B2,F2);
end