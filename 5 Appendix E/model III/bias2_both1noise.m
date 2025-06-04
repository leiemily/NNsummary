clear;
load('bias2-10tpn-step1.mat');
num_step = 1:5892; %前6/10
S_step = S_save2(num_step);
F_step = F_save2(num_step);
dsdt_step = zeros(1,length(num_step));
dFdt_step = dFdt_save2(num_step);
ds2dt2_step = zeros(1,length(num_step));
dF2dt2_step = dF2dt2_save2(num_step);

load('bias2-10tpn-mix1.mat');
num_sin1 = 1:8;%前8/10
s1 = reshape(S_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
F1 = reshape(F_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
dsdt1 = reshape(dsdt_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
dFdt1 = reshape(dFdt_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
ds2dt21 = reshape(ds2dt2_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);
dF2dt21 = reshape(dF2dt2_save(num_sin1,:)',1,length(num_sin1)*discreteT_length);

S_sin1 = s1;
F_sin1 = F1;
dsdt_sin1 = dsdt1;
dFdt_sin1 = dFdt1;
ds2dt2_sin1 = ds2dt21;
dF2dt2_sin1 = dF2dt21;

aS_save2 = [S_step,S_sin1];
aF_save2 = [F_step,F_sin1];
adsdt_save2 = [dsdt_step,dsdt_sin1];
adFdt_save2 = [dFdt_step,dFdt_sin1];
ads2dt2_save2 = [ds2dt2_step,ds2dt2_sin1];
adF2dt2_save2 = [dF2dt2_step,dF2dt2_sin1];

%% s f
S_lin = 0:0.04:2; Ls = length(S_lin);
F_lin = 0:0.02:1; LF = length(F_lin);
[aF1,aS1] = meshgrid(F_lin,S_lin);
%%
save('bias2-10tpn-both1-noise.mat','aS_save2','aF_save2','adsdt_save2','adFdt_save2',...
     'ads2dt2_save2','adF2dt2_save2','aF1','aS1','F_lin','S_lin');
