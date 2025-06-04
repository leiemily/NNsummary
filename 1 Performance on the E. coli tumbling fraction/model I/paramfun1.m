function eq = paramfun1(x,c)
%% parameters
kA = 100;
kY = 130;
kB = 7.5;
kZ = 8.45;

muY = 0.1;
muB = 1;

KB = 0.25;
KBp = 6.5;
KY = 0.65;
KZ = 1;

Bt = 2;
Tt = 5/3;
Yt = 18;
% Zt = 1.1;
Zt = 1.23;
%% x = [B Y Z Tp Bp Yp]
eq(1) = kA * (Tt * c - x(4)) - kY * x(2) * x(4) - kB * x(1) * x(4);
eq(2) = kY * x(2) * x(4) - muY * x(6) -kZ * x(3) * x(6);
eq(3) = kB * x(1) * x(4) - muB * x(5);
eq(4) = x(2) - (Yt - (1 + KZ * x(3)) * x(6)) / (1 + KY * x(4));
eq(5) = x(3) - Zt/(1 + KZ * x(6));
eq(6) = x(1) - (Bt - (1 + KBp * Tt * c) * x(5)) / (1 + KB * x(4));