function sigma = runadsim(Cs,MW,model,params,S)

sigma = zeros(1,length(Cs));
h = waitbar(0,'Running Adsorption Simulation');
for k = 1:length(Cs)
    [coverage, ~] = adsim(MW,Cs(k),model,params,S);
    sigma(k) = mean(coverage)/512^2;
    waitbar(k/length(Cs));
end
close(h)
clear h