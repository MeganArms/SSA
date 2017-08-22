function mu = getonrate(obj)

C = connectEvents(obj,2);
[count, ~] = antiResidenceTimeStat(C,'ExposureTime',1);
[T,SF,SEM] = survivalFunction(length(obj.Frame),1,count);
T = T(SF ~= 0); SF = SF(SF ~= 0); SEM = SEM(SF ~= 0);
fitopt = fitoptions('exp1');
fitopt.Weights = SEM;
[f, ~] = fit(T',SF','exp1',fitopt);
mu = -1/f.b;
end