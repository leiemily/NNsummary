clear;
%% main
j = 0;

%% paths
folderPath_varieddata = '.\dataset1\';
folderPath_NNdata = '.\dataset2\';

trainsetName0 = ['exp1-',num2str(j),'tp-save1d-noise'];
load(fullfile(folderPath_NNdata, trainsetName0));

testsetName0 = ['exp1-',num2str(j),'tp-test1d-noise'];
savePath = fullfile(folderPath_NNdata, testsetName0);

%% test errs
figure(1);
subplot(2,2,1);
plot(ep_number, vali_error1,'o-');
set(gca, 'YScale', 'log');

subplot(2,2,2);
plot(ep_number, test_error1,'o-');
set(gca, 'YScale', 'log');

subplot(2,2,3);
plot(ep_number, vali_error2,'o-');
set(gca, 'YScale', 'log');

subplot(2,2,4);
plot(ep_number, test_error2,'o-');
set(gca, 'YScale', 'log');

loss_errend = [vali_error1(end),vali_error2(end),test_error1(end),test_error2(end)];

%%
datasetName0 = ['exp1-',num2str(j),'tp-both1-noise'];
load(fullfile(folderPath_varieddata, datasetName0),'S_lin','F_lin');
simu1 = Y_surf1;
fit_Gsurfaces;

%%
datasetName1 = ['exp1-',num2str(j),'tp-step1'];
load(fullfile(folderPath_varieddata, datasetName1));

fit1d_step1;

vali_set1 = 21:25;
test_set1 = 26:28;
S_save_step = S_save;
F_pre_step = F_pre1;
F_step = Fsmooth_save2;
t00_step = t00;
Range_step = range_save;

num_test_set1 = Length2(test_set1(1))+1:Length2(test_set1(end)+1);
F_step_test = F_step(num_test_set1);
Fpre_step_test = F_pre_step(num_test_set1);

num_vali_set1 = Length2(vali_set1(1))+1:Length2(vali_set1(end)+1);
F_step_vali = F_step(num_vali_set1);
Fpre_step_vali = F_pre_step(num_vali_set1);

train_set1 = 1:20;
num_train_set1 = Length2(train_set1(1))+1:Length2(train_set1(end)+1);
dFtrain_step = dF_train1(num_train_set1);
dF_step = dFdtsmooth_save2(num_train_set1);

%%
datasetName2 = ['exp1-',num2str(j),'tp-mix1'];
load(fullfile(folderPath_varieddata, datasetName2));
fit1d_mix1;

% vali_set2 = 0;
test_set2 = 7:9;
S_save_mix = S_save;
% dS_save_mix = dsdt_save/t_scale;
F_pre_mix = F_pre1;
F_mix = Fsmooth_save;
Fnoise_mix = Fnoise_save;
t00_mix = t00;
Range_mix = range_save;

F_mix_test = reshape(F_mix(test_set2,:)',[1,length(test_set2)*discreteT_length]);
Fpre_mix_test = reshape(F_pre_mix(test_set2,:)',[1,length(test_set2)*discreteT_length]);

% F_mix_vali = reshape(F_mix(vali_set2,:)',[1,length(vali_set2)*discreteT_length]);
% Fpre_mix_vali = reshape(F_pre_mix(vali_set2,:)',[1,length(vali_set2)*discreteT_length]);

train_set2 = 1:6;
dF_mix = reshape(dFdtsmooth_save(train_set2,:)',[1,length(train_set2)*discreteT_length]);
dFtrain_mix = reshape(dF_train1(train_set2,:)',[1,length(train_set2)*discreteT_length]);


%% calculate err
%% train
N_train_step = length(dF_step);
N_train_mix = length(dF_mix);
train_all = [dF_step,dF_mix];
trainpre_all = [dFtrain_step,dFtrain_mix];
N_train_all = length(train_all);

error_RAEtrain_step = sum(abs(dF_step- dFtrain_step))/sum(abs(dF_step));
error_RSEtrain_step = sum((dF_step- dFtrain_step).^2)/sum((dF_step).^2);

error_RAEtrain_mix = sum(abs(dF_mix- dFtrain_mix))/sum(abs(dF_mix));
error_RSEtrain_mix = sum((dF_mix- dFtrain_mix).^2)/sum((dF_mix).^2);

error_RAEtrain_all = sum(abs(train_all- trainpre_all))/sum(abs(train_all));
error_RSEtrain_all = sum((train_all- trainpre_all).^2)/sum((train_all).^2);

%% vali
N_vali_step = length(F_step_vali);
N_vali_mix = 0;

vali_all = F_step_vali;
valipre_all = Fpre_step_vali;
N_vali_all = length(vali_all);

error_MAEvali_step = 1/N_vali_step * sum(abs(F_step_vali- Fpre_step_vali)./abs(F_step_vali));
error_MSEvali_step = 1/N_vali_step * sum((F_step_vali- Fpre_step_vali).^2./(F_step_vali).^2);
error_R2vali_step = 1- sum((F_step_vali- Fpre_step_vali).^2)/sum((F_step_vali-mean(F_step_vali)).^2);

error_MAEvali_mix = 0;
error_MSEvali_mix = 0;
error_R2vali_mix = 0;

error_MAEvali_all = 1/N_vali_all * sum(abs(vali_all- valipre_all)./abs(vali_all));
error_MSEvali_all = 1/N_vali_all * sum((vali_all- valipre_all).^2./(vali_all).^2);
error_R2vali_all = 1- sum((vali_all- valipre_all).^2)/sum((vali_all-mean(vali_all)).^2);

%% test
N_test_step = length(F_step_test);
N_test_mix = length(F_mix_test);

test_all = [F_step_test,F_mix_test];
testpre_all = [Fpre_step_test,Fpre_mix_test];
N_test_all = length(test_all);

error_MAEtest_step = 1/N_test_step * sum(abs(F_step_test- Fpre_step_test)./abs(F_step_test));
error_MSEtest_step = 1/N_test_step * sum((F_step_test- Fpre_step_test).^2./(F_step_test).^2);
error_R2test_step = 1- sum((F_step_test- Fpre_step_test).^2)/sum((F_step_test-mean(F_step_test)).^2);

error_MAEtest_mix = 1/N_test_mix * sum(abs(F_mix_test- Fpre_mix_test)./abs(F_mix_test));
error_MSEtest_mix = 1/N_test_mix * sum((F_mix_test- Fpre_mix_test).^2./(F_mix_test).^2);
error_R2test_mix = 1- sum((F_mix_test- Fpre_mix_test).^2)/sum((F_mix_test-mean(F_mix_test)).^2);

error_MAEtest_all = 1/N_test_all * sum(abs(test_all- testpre_all)./abs(test_all));
error_MSEtest_all = 1/N_test_all * sum((test_all- testpre_all).^2./(test_all).^2);
error_R2test_all = 1- sum((test_all- testpre_all).^2)/sum((test_all-mean(test_all)).^2);

%%
err_all = [[N_train_step,error_RAEtrain_step,error_RSEtrain_step,0];
    [N_train_mix,error_RAEtrain_mix,error_RSEtrain_mix,0];
    [N_train_all,error_RAEtrain_all,error_RSEtrain_all,0];
    [N_vali_step,error_MAEvali_step,error_MSEvali_step,error_R2vali_step];
    [N_vali_mix,error_MAEvali_mix,error_MSEvali_mix,error_R2vali_mix];
    [N_vali_all,error_MAEvali_all,error_MSEvali_all,error_R2vali_all];
    [N_test_step,error_MAEtest_step,error_MSEtest_step,error_R2test_step];
    [N_test_mix,error_MAEtest_mix,error_MSEtest_mix,error_R2test_mix];
    [N_test_all,error_MAEtest_all,error_MSEtest_all,error_R2test_all]]';
 

%% save
save(savePath,'dt00','t_scale','err_all','ep_number','loss_errend');
save(savePath,'test_set1','test_error1','vali_error1','err1','t00_step','T_save2','Length2',...
    'F_step','F_pre_step','S_save_step','Range_step','-append');
save(savePath,'test_set2','test_error2','vali_error2','err2','t00_mix',...
    'Fnoise_mix','F_mix','F_pre_mix','S_save_mix','Range_mix','-append');
