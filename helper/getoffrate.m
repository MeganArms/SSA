function k = getoffrate(obj)

C = getCoordinates(obj,'yes');
[count, edges] = ResidenceTimeStat(C,'ExposureTime',1);
[T,SF,SEM] = survivalFunction(length(obj.Frame),1,count);
T = T(SF ~= 0); SF = SF(SF ~= 0); SEM = SEM(SF ~= 0);
fitopt = fitoptions('exp1');
fitopt.Weights = SEM;
[f, gof] = fit(T',SF','exp1',fitopt);
k = -1/f.b;
end