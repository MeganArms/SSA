function viewtraces(traces,Nframes,et)
figure;
t = et:et:Nframes*et;
for i = 1:size(traces,1)
    plot(t,traces(i,:));title(['Trace ',num2str(i)]);
    pause(1)
end