function b_i = mlesolver(lag,Nframes)

syms t c1 l1 c2 l2 l3
b = [c1 c2 l1 l2 l3];
f(t) = c1*l1*exp(-l1*t) + c2*l2*exp(-l2*t) + (1-c1-c2)*exp(-l3*t);
logf(t) = log(f); % Then I'm supposed to sum over the observations
LOGF = matlabFunction(log(f),'Vars',{t, b});
tt = (lag:lag:lag*Nframes);

% Uc1 = diff(logf,c1); Ul1 = diff(logf,l1); Uc2 = diff(logf,c2); Ul2 = diff(logf,l2);
% Ul3 = diff(logf,l3);
% Ic1c1 = diff(Uc1,c1); Ic1l1 = diff(Uc1,l1); Ic1c2 = diff(U
Usym = [diff(logf,c1),diff(logf,l1),diff(logf,c2),diff(logf,l2),diff(logf,l3)];
I1 = diff(Usym,c1); I2 = diff(Usym,l1); I3 = diff(Usym,c2); I4 = diff(Usym,l2);
I5 = diff(Usym,l3);
Isym = [I1,I2,I3,I4,I5];

UFunc = matlabFunction(Usym,'Vars',{t b});
IFunc = matlabFunction(Isym,'Vars',{t b});

nparams = length(b);
b_i = [0.5 0.1 0.5 0.01 1];
nll_next = -sum(LOGF(tt,b_i));
nIter = 1; iterMax = 10^6; N = length(tt);
while nIter < iterMax
    nll_prev = nll_next;
    U = reshape(sum(reshape(UFunc(tt,b_i),[N,1,nparams])),[nparams,1]);
    I = reshape(mean(reshape(IFunc(tt,b_i),[N,nparams,nparams])),[nparams, nparams])';
    b_current = b_i + (I\U)';
    nll_next = -2*sum(LOGF(tt,b_current));
    if nll_next > nll_prev
        b_i = (b_i + b_current)./2;
        nIter = nIter + 1;
    elseif nll_next < nll_prev
        b_i = b_current;
        nIter = nIter + 1;    
    elseif (nll_prev - nll_next) < 0.025
        break
    elseif nIter == iterMax
        disp('Max iterations reached');
        return
    else
        disp('something else is happening');
    end
end