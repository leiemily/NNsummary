clear;
%% main
j = 0;

%% paths
folderPath_varieddata = '.\dataset1\';
folderPath_NNdata = '.\dataset2\';

trainsetName0 = ['exp1-',num2str(j),'tp-save2d-noise'];
load(fullfile(folderPath_NNdata, trainsetName0));

testsetName0 = ['exp1-',num2str(j),'tp-test2d-noise'];
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
datasetName1 = ['exp1-',num2str(j),'tp-step1'];
load(fullfile(folderPath_varieddata, datasetName1));
fit2d_trapz1;

vali_set1 = 21:25;
test_set1 = 26:28;
S_save_trapz = S_save;
F_pre_trapz = F_pre1;
F_trapz = Fsmooth_save2;
t00_trapz = t00;
Range_trapz = range_save;

num_test_set1 = Length2(test_set1(1))+1:Length2(test_set1(end)+1);
F_trapz_test = F_trapz(num_test_set1);
Fpre_trapz_test = F_pre_trapz(num_test_set1);

num_vali_set1 = Length2(vali_set1(1))+1:Length2(vali_set1(end)+1);
F_trapz_vali = F_trapz(num_vali_set1);
Fpre_trapz_vali = F_pre_trapz(num_vali_set1);

train_set1 = 1:20;
num_train_set1 = Length2(train_set1(1))+1:Length2(train_set1(end)+1);
dF2train_trapz = dF2_train1(num_train_set1);
dF2_trapz = dF2dt2smooth_save2(num_train_set1);

%%
datasetName2 = ['exp1-',num2str(j),'tp-mix1'];
load(fullfile(folderPath_varieddata, datasetName2));
fit2d_mix1;

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
dF2train_mix = reshape(dF2_train1(train_set2,:)',[1,length(train_set2)*discreteT_length]);
dF2_mix = reshape(dF2dt2smooth_save(train_set2,:)',[1,length(train_set2)*discreteT_length]);


%% calculate err
%% train
N_train_trapz = length(dF2_trapz);
N_train_mix = length(dF2_mix);
train_all = [dF2_trapz,dF2_mix];
trainpre_all = [dF2train_trapz,dF2train_mix];
N_train_all = length(train_all);

error_RAEtrain_trapz = sum(abs(dF2_trapz- dF2train_trapz))/sum(abs(dF2_trapz));
error_RSEtrain_trapz = sum((dF2_trapz- dF2train_trapz).^2)/sum((dF2_trapz).^2);

error_RAEtrain_mix = sum(abs(dF2_mix- dF2train_mix))/sum(abs(dF2_mix));
error_RSEtrain_mix = sum((dF2_mix- dF2train_mix).^2)/sum((dF2_mix).^2);

error_RAEtrain_all = sum(abs(train_all- trainpre_all))/sum(abs(train_all));
error_RSEtrain_all = sum((train_all- trainpre_all).^2)/sum((train_all).^2);

%% vali
N_vali_trapz = length(F_trapz_vali);
N_vali_mix = 0;

vali_all = F_trapz_vali;
valipre_all = Fpre_trapz_vali;
N_vali_all = length(vali_all);

error_MAEvali_trapz = 1/N_vali_trapz * sum(abs(F_trapz_vali- Fpre_trapz_vali)./abs(F_trapz_vali));
error_MSEvali_trapz = 1/N_vali_trapz * sum((F_trapz_vali- Fpre_trapz_vali).^2./(F_trapz_vali).^2);
error_R2vali_trapz = 1- sum((F_trapz_vali- Fpre_trapz_vali).^2)/sum((F_trapz_vali-mean(F_trapz_vali)).^2);

error_MAEvali_mix = 0;
error_MSEvali_mix = 0;
error_R2vali_mix = 0;

error_MAEvali_all = 1/N_vali_all * sum(abs(vali_all- valipre_all)./abs(vali_all));
error_MSEvali_all = 1/N_vali_all * sum((vali_all- valipre_all).^2./(vali_all).^2);
error_R2vali_all = 1- sum((vali_all- valipre_all).^2)/sum((vali_all-mean(vali_all)).^2);

%% test
N_test_trapz = length(F_trapz_test);
N_test_mix = length(F_mix_test);

test_all = [F_trapz_test,F_mix_test];
testpre_all = [Fpre_trapz_test,Fpre_mix_test];
N_test_all = length(test_all);

error_MAEtest_trapz = 1/N_test_trapz * sum(abs(F_trapz_test- Fpre_trapz_test)./abs(F_trapz_test));
error_MSEtest_trapz = 1/N_test_trapz * sum((F_trapz_test- Fpre_trapz_test).^2./(F_trapz_test).^2);
error_R2test_trapz = 1- sum((F_trapz_test- Fpre_trapz_test).^2)/sum((F_trapz_test-mean(F_trapz_test)).^2);

error_MAEtest_mix = 1/N_test_mix * sum(abs(F_mix_test- Fpre_mix_test)./abs(F_mix_test));
error_MSEtest_mix = 1/N_test_mix * sum((F_mix_test- Fpre_mix_test).^2./(F_mix_test).^2);
error_R2test_mix = 1- sum((F_mix_test- Fpre_mix_test).^2)/sum((F_mix_test-mean(F_mix_test)).^2);

error_MAEtest_all = 1/N_test_all * sum(abs(test_all- testpre_all)./abs(test_all));
error_MSEtest_all = 1/N_test_all * sum((test_all- testpre_all).^2./(test_all).^2);
error_R2test_all = 1- sum((test_all- testpre_all).^2)/sum((test_all-mean(test_all)).^2);

%%
err_all = [[N_train_trapz,error_RAEtrain_trapz,error_RSEtrain_trapz,0];
    [N_train_mix,error_RAEtrain_mix,error_RSEtrain_mix,0];
    [N_train_all,error_RAEtrain_all,error_RSEtrain_all,0];
    [N_vali_trapz,error_MAEvali_trapz,error_MSEvali_trapz,error_R2vali_trapz];
    [N_vali_mix,error_MAEvali_mix,error_MSEvali_mix,error_R2vali_mix];
    [N_vali_all,error_MAEvali_all,error_MSEvali_all,error_R2vali_all];
    [N_test_trapz,error_MAEtest_trapz,error_MSEtest_trapz,error_R2test_trapz];
    [N_test_mix,error_MAEtest_mix,error_MSEtest_mix,error_R2test_mix];
    [N_test_all,error_MAEtest_all,error_MSEtest_all,error_R2test_all]]';
 
%% save
save(savePath,'dt00','t_scale','err_all','ep_number','loss_errend');
save(savePath,'test_set1','test_error1','vali_error1','err1','t00_trapz','T_save2','Length2',...
    'F_trapz','F_pre_trapz','S_save_trapz','Range_trapz','-append');
save(savePath,'test_set2','test_error2','vali_error2','err2','t00_mix',...
    'Fnoise_mix','F_mix','F_pre_mix','S_save_mix','Range_mix','-append');