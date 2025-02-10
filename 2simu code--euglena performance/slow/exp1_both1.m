clear;
%%
load('exp1-0tp-step1.mat');
num_step = 1:7828;% 前20/28
S_step = S_save2(num_step);
F_step = Fsmooth_save2(num_step);
dsdt_step = zeros(1,length(num_step));
dFdt_step = dFdtsmooth_save2(num_step);
ds2dt2_step = zeros(1,length(num_step));
dF2dt2_step = dF2dt2smooth_save2(num_step);

%%
load('exp1-0tp-mix1.mat');
num_sin1 = 1:6;%前6/9
s1 = reshape(S_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
F1 = reshape(Fsmooth_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
dsdt1 = reshape(dsdt_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
dFdt1 = reshape(dFdtsmooth_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
ds2dt21 = reshape(ds2dt2_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
dF2dt21 = reshape(dF2dt2smooth_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);

S_sin1 = s1;
F_sin1 = F1;
dsdt_sin1 = dsdt1;
dFdt_sin1 = dFdt1;
ds2dt2_sin1 = ds2dt21;
dF2dt2_sin1 = dF2dt21;
%%
aS_save2 = [S_step,S_sin1];
aF_save2 = [F_step,F_sin1];
adsdt_save2 = [dsdt_step,dsdt_sin1];
adFdt_save2 = [dFdt_step,dFdt_sin1];
ads2dt2_save2 = [ds2dt2_step,ds2dt2_sin1];
adF2dt2_save2 = [dF2dt2_step,dF2dt2_sin1];

%% s f
S_lin = 0:0.01:0.5; Ls = length(S_lin);
F_lin = 0:0.01:0.3; LF = length(F_lin);
[aF1,aS1] = meshgrid(F_lin,S_lin);

%%
save('exp1-0tp-both1-noise.mat','aS_save2','aF_save2','adsdt_save2','adFdt_save2',...
     'ads2dt2_save2','adF2dt2_save2','aF1','aS1','F_lin','S_lin');