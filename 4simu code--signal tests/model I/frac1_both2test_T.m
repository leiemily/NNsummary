clear;

load('frac1-1tp-trapz1.mat');
num_trapz1 = 1:14014;
S_trapz1 = S_save2(num_trapz1);
F_trapz1 = F_save2(num_trapz1);
dsdt_trapz1 = dsdt_save2(num_trapz1);
dFdt_trapz1 = dFdt_save2(num_trapz1);
ds2dt2_trapz1 = zeros(1,length(num_trapz1));
dF2dt2_trapz1 = dF2dt2_save2(num_trapz1);

%%
aS_save2 = S_trapz1;
aF_save2 = F_trapz1;
adsdt_save2 = dsdt_trapz1;
adFdt_save2 = dFdt_trapz1;
ads2dt2_save2 = ds2dt2_trapz1;
adF2dt2_save2 = dF2dt2_trapz1;
%% s f
S_lin = 0:0.01:0.5; Ls = length(S_lin);
F_lin = 0:0.01:0.5; LF = length(F_lin);
[aF1,aS1] = meshgrid(F_lin,S_lin);

%%
save('frac1-1tp-both1-freeT.mat','aS_save2','aF_save2','adsdt_save2','adFdt_save2',...
     'ads2dt2_save2','adF2dt2_save2','aF1','aS1');