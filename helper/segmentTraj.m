function [trajs, ctrajs] = segmentTraj(spotEvents)

trajs = struct;
ctrajs = struct;
l = 1;
m = 1;
for k = 1:length(spotEvents)
    if spotEvents(k).eventfreq == 1
        trajs(l).trajectory = spotEvents(k).trajectory;
        trajs(l).coordinates = spotEvents(k).coordinates;
        trajs(l).brightness = spotEvents(k).brightness;
        trajs(l).frames = spotEvents(k).frames;
        trajs(l).eventfreq = spotEvents(k).eventfreq;
        trajs(l).std = spotEvents(k).std;
        l = l + 1;
    else
        ctrajs(m).trajectory = spotEvents(k).trajectory;
        ctrajs(m).coordinates = spotEvents(k).coordinates;
        ctrajs(m).brightness = spotEvents(k).brightness;
        ctrajs(m).frames = spotEvents(k).frames;
        ctrajs(m).eventfreq = spotEvents(k).eventfreq;
        ctrajs(m).std = spotEvents(k).std;
        m = m + 1;
    end
end