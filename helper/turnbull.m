% Turnbull Intervals
a = 0;
b = 60*60*24;
L = cellfun(@(x)x(1),F);
R = cellfun(@(x)x(end),F);
cens_int = L == 1;
cens_r = R == Nframes*lag;
deadtime = t1-a; % Time between experiment start and obs start
L(cens_int) = L(cens_int)+deadtime;
R(cens_r) = Inf;

t = 0:2:Nframes*lag; t = [t,Inf];
p = zeros(length(L),length(t));
for k = 1:length(L)
    Li = find(t==L(k));
    Ri = find(t==R(k));
    w = 1/(Ri-Li);
    p(k,Li:Ri) = repmat(w,length(Li:Ri),1);
end

