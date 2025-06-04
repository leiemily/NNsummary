function eq = paramfun_m(x,c)
%% parameters
kR = 3.82e-2;
a0 = 0.5;
m0 = 1;
Nr = 6;
alpha0 = 1.7;
KI = 18.2;
KA = 3000;


%% 
eq = kR * (1 - (1 + exp (Nr * (- alpha0 *(x - m0) + log((1 + c/KI)/(1 + c/KA))))).^(-1)/a0);