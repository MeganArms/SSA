function trace = picktraces(trkid,Nframes,F1,B)

trace = zeros(Nframes,1);
k = 1; trknum = 1;
while k <= Nframes && trknum <= length(trkid)
    b = B{trkid(trknum)};
    f = F1{trkid(trknum)};
    if k == f(1)
        trace(f) = b;
        k = k+length(f);
        trknum = trknum + 1;   
    else
        inds = k:f(1)-1;
        trace(inds) = zeros(length(inds),1);
        k = f(1);
    end
end
if k < Nframes
    trace(k:Nframes) = zeros(length(k:Nframes),1);
end