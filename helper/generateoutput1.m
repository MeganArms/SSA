function [survFuncs,intensities,numbers] = generateoutput1(et,Nframes,C1,F1,spotEvents1,objs_link1,eventfreq1,varargin)

[sfc,tc,semc,bc,numec] = crtdPB(spotEvents1,Nframes,et,objs_link1);
[icounts,dtime,firstframe,lastframe,~,avgb,rangeb,maxb,rcdfb] = factorsPB(spotEvents1,Nframes,objs_link1);
[f,x,~,~,~,CC,b,xx,nume] = residenceTime(C1,F1,Nframes,et);
[fc,xc,reside,cens] = kmPB(spotEvents1,Nframes,et,objs_link1);
num1 = sum(firstframe == 1); 
fffract = num1/length(firstframe);
numlast = sum(lastframe == Nframes);
lffract = numlast/length(lastframe);
numtrajc = length(spotEvents1);
numobj = size(objs_link1,2);
numtraj = length(C1);
xc(1) = 0; x(1) = 0;

if isempty(varargin)
    survFuncs = {sfc,tc,semc,bc,numec,fc,xc,reside,cens,CC,xx,b,nume,f,x};
    intensities = {icounts(:,1),icounts(:,2),icounts(:,3),icounts(:,4),dtime,avgb,rangeb,maxb,rcdfb};
    numbers = {firstframe,lastframe,num1,fffract,numlast,lffract,numtrajc,numtraj,numobj,eventfreq1,Nframes};
else
    survFuncs = varargin{1};
    newIndex = size(survFuncs,1)+1;
    intensities = varargin{2};
    numbers = varargin{3};
    
    survFuncs(newIndex,:) = {sfc,tc,semc,bc,numec,fc,xc,reside,cens,CC,xx,b,nume,f,x};
    intensities(newIndex,:) = {icounts(:,1),icounts(:,2),icounts(:,3),icounts(:,4),dtime,avgb,rangeb,maxb,rcdfb};
    numbers(newIndex,:) = {firstframe,lastframe,num1,fffract,numlast,lffract,numtrajc,numtraj,numobj,eventfreq1,Nframes};
end

% objs_link1 = objs_link;
% clear objs_link objs
% [C1,B1,F1,M1,S1,R1,Frame1] = collate(objs_link1);
% [spotEvents1, objtrajinf, eventfreq1, stuckfreq, empties] = eventlinks2(C1,R1,objs_link1);
% 
% objs_link2 = objs_link;
% clear objs objs_link
% [C2,B2,F2,M2,S2,R2,Frame2] = collate(objs_link2);
% [spotEvents2, objtrajinf, eventfreq2, stuckfreq, empties] = eventlinks2(C2,R2,objs_link2);
end