function eq = paramfun_m(x,c)
%% parameters
m0 = 1;
KI = 18.2;
KA = 3000;

Nr = 6;
alpha0 = 1.7;


V_R = 0.01;
V_B = 0.02;
V_on = 0.02;
V_off = 0.0194;
FliM_max = 45;


%% 
eq(1) = V_R*(1 - (1 + exp (Nr * (- alpha0 *(x(1) - m0) + log((1 + c(1)/KI)/(1 + c(1)/KA))))).^(-1)) -...
    V_B*(1 + exp (Nr * (- alpha0 *(x(1) - m0) + log((1 + c(1)/KI)/(1 + c(1)/KA))))).^(-1);
eq(2) = V_on * (FliM_max - x(2)) - V_off * (1 + exp (Nr * (- alpha0 *(x(1) - m0) + log((1 + c(1)/KI)/(1 + c(1)/KA))))).^(-1)* x(2);