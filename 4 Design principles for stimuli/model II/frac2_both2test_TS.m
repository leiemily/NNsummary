clear;
load('frac2-1tp-step1.mat');
num_step = 1:7007; 
S_step = S_save2(num_step);
F_step = F_save2(num_step);
dsdt_step = dsdt_save2(num_step);
dFdt_step = dFdt_save2(num_step);
ds2dt2_step = zeros(1,length(num_step));
dF2dt2_step = dF2dt2_save2(num_step);

load('frac2-1tp-trapz1.mat');
num_trapz1 = 1:7007;
S_trapz1 = S_save2(num_trapz1);
F_trapz1 = F_save2(num_trapz1);
dsdt_trapz1 = dsdt_save2(num_trapz1);
dFdt_trapz1 = dFdt_save2(num_trapz1);
ds2dt2_trapz1 = zeros(1,length(num_trapz1));
dF2dt2_trapz1 = dF2dt2_save2(num_trapz1);

%%
aS_save2 = [S_step,S_trapz1];
aF_save2 = [F_step,F_trapz1];
adsdt_save2 = [dsdt_step,dsdt_trapz1];
adFdt_save2 = [dFdt_step,dFdt_trapz1];
ads2dt2_save2 = [ds2dt2_step,ds2dt2_trapz1];
adF2dt2_save2 = [dF2dt2_step,dF2dt2_trapz1];

%% s f
S_lin = 0:0.04:2; Ls = length(S_lin);
F_lin = 0:0.01:0.5; LF = length(F_lin);
[aF1,aS1] = meshgrid(F_lin,S_lin);
%%
save('frac2-1tp-both1-freeTS.mat','aS_save2','aF_save2','adsdt_save2','adFdt_save2',...
     'ads2dt2_save2','adF2dt2_save2','aF1','aS1','F_lin','S_lin');