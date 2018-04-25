function [rpb, cpb] = repeatltspb(spotEvents,eventfreq,Nframes)
% Find the lifetimes of the individual events within multiple events

rpb = zeros(length(spotEvents)*max(eventfreq),1);
l = 1;
for k = 1:length(spotEvents)
    lts = regionprops(~isnan(spotEvents(k).trajectory),'Area');
    ltsvec = [lts.Area]';
    if isnan(spotEvents(k).trajectory(1))
        continue
    end
    censtmp = zeros(length(ltsvec),1);
    if spotEvents(k).trajectory(end) == Nframes
        censtmp(end) = 1;
    end
    rpb(l:length(ltsvec)+l-1) = ltsvec;
    cpb(l:length(ltsvec)+l-1) = censtmp;
    l = l + length(lts);
end
rpb = rpb(rpb>0);
cpb = cpb(rpb>0);