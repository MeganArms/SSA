function [coverage,state] = adsim(MW,C,model,params,S)

% Molecular weight in g/mol
% Concentration in mol/L

a = 1;
b = 512*512;
A_FOV = 81.92*81.92/(10^6)^2; % m^2
V_FC = 8*22*.1/(10^3)^3; % m^3
N = 6.022E23; % Avogardro's number
R = 8.314; % J/K/mol
m_p = MW/1000; % kg/mol
Nframes = 3601; % s
% S = 5E-4; % Sticking probability from ref. 20 in Katira 2009 - chose 10x larger
Nb = round(C*V_FC/10^-3*N); % Number of molecules in bulk initially

state = zeros(b,Nframes); % Location changes over i. Time changes over j.
coverage = zeros(Nframes,1);

for k = 1:Nframes
    % Calculate the number of surface molecules with current bulk
    % concentration and temperature
    T = 298 + normrnd(0,0.5); % K
    J = C*sqrt(R*T/m_p); % mol/s/m^2
    Ns = round(S*J*N*A_FOV);
    % Choose random locations for spots
    spots = round((b-a).*rand(Ns,1)+a);
    % Do not allow more molecules to adsorb than were present in bulk
    if Ns > Nb
        Ns = Nb;
    end
    % Generate residence times
    if strcmp(model,'exp1')
        RTs = exprnd(params,Ns,1);
    elseif strcmp(model,'gp')
        RTs = gprnd(params(1),params(2),params(3),[1,Ns]); % k, sigma, theta
    end
    % Find integer time of desorption
    des = k + round(RTs);
    % Do not allow for times longer than the video
    if sum(des > Nframes) > 0
        des(des > Nframes) = Nframes;
    end
    % Do not allow multiple events in the same spot -- slow
%     currstate = state(:,k);
%     overlap = ismember(spots,find(currstate));
%     numoverlap = sum(overlap);
%     spots = spots(~overlap);
%     Ns = Ns - numoverlap;
%     Ndes = numoverlap;
    % Add a tick to locations (i) at time (j) 
    Ndes = 0;
    for l = 1:Ns
        if state(spots(l),k) ~= 1
            state(spots(l),k:des(l)) = state(spots(l),k:des(l)) + 1;
        else
            Ndes = Ndes + 1;
        end
    end
%     state(spots,k:des) = state(spots,k:des) + 1;
    coverage(k) = sum(state(:,k) > 0);
    % Find the number desorbing at k
    if k == 1
        Ndes = 0;
    else
        Ndes = Ndes + sum(state(:,k)-state(:,k-1) < 0);
    end
    % Calculate the new bulk number of molecules
    Nb = C*V_FC*N - Ns + Ndes;
    % Calculate the concentration after change in Nb
    C = Nb/N/V_FC; % mol/m^3
end
    
    