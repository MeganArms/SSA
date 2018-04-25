function [alpha, C] = plmle(r,c)

uncens = r(~c);
alpha = sum(uncens)/sum(log(r));
C = (alpha - 1)*2^(alpha-1);