%% fit surfaces
Ls = length(S_lin);
LF = length(F_lin);
[F1,S1] = meshgrid(F_lin,S_lin);
%% G1 pred
G1_nn = reshape(simu1(:,1), LF,Ls)';

%% G2 pred
G2_nn = reshape(simu1(:,2), LF,Ls)';

%% figures
figure(2);
subplot(1,2,1);
pcolor(F1,S1,G1_nn);
shading interp; 
colorbar; colormap(jet); 
xlabel('$f$','Interpreter','latex');
ylabel('$s$','Interpreter','latex');
title('Learned $G_1(f,s)$','Interpreter','latex');

subplot(1,2,2);
pcolor(F1,S1,G2_nn);
shading interp; 
colorbar; colormap(jet); 
xlabel('$f$','Interpreter','latex');
ylabel('$s$','Interpreter','latex');
title('Learned $G_2(f,s)$','Interpreter','latex');