function F = objectloss(s,Nframes,et)

tmax = (Nframes - 2)*et; T = et:et:tmax;
nmax = 2*(Nframes-2);
% timevec = 0:et:2*tmax; inds = 1:length(timevec);
% p = 0.029/4.3*exp(-timevec/4.3)+0.17/0.56*exp(-timevec/0.56)+0.45/0.08*exp(-timevec/0.08);

timevec = et:et:2*tmax;
p = 0.029/4.3*exp(-timevec/4.3)+0.17/0.56*exp(-timevec/0.56)+0.45/0.08*exp(-timevec/0.08);
% p = @(t)0.029/4.3*exp(-t/4.3)+0.17/0.56*exp(-t/0.56)+0.45/0.08*exp(-t/0.08);
q = zeros(Nframes-1,1);
ns = tril(meshgrid(1:nmax)'); % changed 2 to 1
gs = meshgrid(0:nmax-1); % changed 2:nmax to 0:nmax-2
r = 1:nmax; % r = [1,1:(nmax-1)];
nc = tril(squareform(pdist(r')));
nc = [nc(2:end,:);nmax*ones(1,nmax)];
% nc = [ones(1,nmax);nc(1:end-1,:)];
nc(nc<1) = 1;

h = waitbar(0,'Analyzing...');
for x = 1:Nframes
    tic;gg = tril(gs,1-x);toc
    tic;nn = tril([zeros(x-1,nmax);ns(x:end,:)],1-x);toc
    tic;bino = (1-s).^(nn-gg).*s.^gg.*(gg+1); toc % LONG
    tic;gg1 = gg; gg1(gg1<2) = 2;toc
    tic;g1 = tril(cumprod(gg1-1,2));toc
    tic;g1(g1<1) = 1;toc
    tic;ncx = tril(cumprod([ones(nmax,x+1),nc(:,x+2:end)],2));toc
    % ncx(ncx<1) = 1;
    tic;nck = ncx./g1;toc
    tic;inds = nn(:,x) > 0;toc
    tic;q(x) = sum(p'.*inds.*sum(bino.*nck,2),1,'omitnan');toc
    waitbar(x/Nframes);
end
close(h); clear h
F = zeros(length(T),1); N = sum(q(2:tmax/et));
for k = 1:length(T)
    t = T(k);
    F(k) = sum(q((1+round(t/et)):round(tmax/et)))/N;
end


