clear;
%% main
j = 1;
plot_number = 1;%plot when ==1 else do not plot
%% paths
folderPath_varieddata = '.\dataset1\';
folderPath_NNdata = '.\dataset2\';

trainsetName0 = ['frac2-',num2str(j),'tp-save2d-freeTS'];
load(fullfile(folderPath_NNdata, trainsetName0));
testsetName0 = ['frac2-',num2str(j),'tp-test2d-freeTS'];
savePath = fullfile(folderPath_NNdata, testsetName0);

% trainsetName0 = ['frac2-',num2str(j),'tp-save2d-freeTSL'];
% load(fullfile(folderPath_NNdata, trainsetName0));
% testsetName0 = ['frac2-',num2str(j),'tp-test2d-freeTSL'];
% savePath = fullfile(folderPath_NNdata, testsetName0);

% trainsetName0 = ['frac2-',num2str(j),'tp-save2d-freeLT'];
% load(fullfile(folderPath_NNdata, trainsetName0));
% testsetName0 = ['frac2-',num2str(j),'tp-test2d-freeLT'];
% savePath = fullfile(folderPath_NNdata, testsetName0);

% trainsetName0 = ['frac2-',num2str(j),'tp-save2d-freeLS'];
% load(fullfile(folderPath_NNdata, trainsetName0));
% testsetName0 = ['frac2-',num2str(j),'tp-test2d-freeLS'];
% savePath = fullfile(folderPath_NNdata, testsetName0);

% trainsetName0 = ['frac2-',num2str(j),'tp-save2d-freeL'];
% load(fullfile(folderPath_NNdata, trainsetName0));
% testsetName0 = ['frac2-',num2str(j),'tp-test2d-freeL'];
% savePath = fullfile(folderPath_NNdata, testsetName0);
%% test errs
% Fi = 2;
% Fj = 2;
% figure(1);
% subplot(Fi,Fj,1);
% plot(ep_number, test_error1,'o-');
% set(gca, 'YScale', 'log');
% 
% subplot(Fi,Fj,2);
% plot(ep_number, test_error2,'o-');
% set(gca, 'YScale', 'log');
% 
% subplot(Fi,Fj,3);
% plot(ep_number, test_error3,'o-');
% set(gca, 'YScale', 'log');

% subplot(Fi,Fj,4);
% plot(ep_number, test_error4,'o-');
% set(gca, 'YScale', 'log');

%%
datasetName1 = ['frac2-',num2str(j),'tp-step1'];
load(fullfile(folderPath_varieddata, datasetName1));
simu_step = Y_all1;
repeat_number1 = 1;
[err1,F_pre1,dF2_pre1,dF2_train1] = fit2d_disconti1(simu_step,Length2,Numberof_F,t_scale,T_save2,...
    F_save2,dFdt_save2,dF2dt2_save2,ds2dt2_save2,Numberof_F/4,4,repeat_number1,plot_number);

train_set1 = 1:6;
test_set1 = 7:16;

S_save_step = S_save;
F_pre_step = F_pre1;
Ffree_step = F_save;
F_step = F_save2;
t00_step = t00;
Range_step = range_save0;

num_test_set1 = Length2(test_set1(1))+1:Length2(test_set1(end)+1);
F_step_test = F_step(num_test_set1);
Fpre_step_test = F_pre_step(num_test_set1);

num_train_set1 = Length2(train_set1(1))+1:Length2(train_set1(end)+1);
dF2train_step = dF2_train1(num_train_set1);
dF2_step = dF2dt2_save2(num_train_set1)/t_scale^2;

%%
datasetName2 = ['frac2-',num2str(j),'tp-trapz1'];
load(fullfile(folderPath_varieddata, datasetName2));
simu_trapz = Y_all2;
repeat_number2 = 2;
[err2,F_pre2,dF2_pre2,dF2_train2] = fit2d_disconti1(simu_trapz,Length2,Numberof_F,t_scale,T_save2,...
    F_save2,dFdt_save2,dF2dt2_save2,ds2dt2_save2,Numberof_F/4,4,repeat_number2,plot_number);

train_set2 = 1:6;
test_set2 = 7:16;

S_save_trapz = S_save;
F_pre_trapz = F_pre2;
Ffree_trapz = F_save;
F_trapz = F_save2;
t00_trapz = t00;
Range_trapz = range_save0;

num_test_set2 = Length2(test_set2(1))+1:Length2(test_set2(end)+1);
F_trapz_test = F_trapz(num_test_set2);
Fpre_trapz_test = F_pre_trapz(num_test_set2);

num_train_set2 = Length2(train_set2(1))+1:Length2(train_set2(end)+1);
dF2train_trapz = dF2_train2(num_train_set2);
dF2_trapz = dF2dt2_save2(num_train_set2)/t_scale^2;


%%
datasetName3 = ['frac2-',num2str(j),'tp-mix1'];
load(fullfile(folderPath_varieddata, datasetName3));
simu_mix = Y_all3;
repeat_number3 = 3;

[err3,F_pre3,dF2_pre3,second_term_test3,dF2_train3] = fit2d_conti1(simu_mix,discreteT_length,Numberof_F,t00,t_scale,...
    F_save,dFdt_save,dF2dt2_save,ds2dt2_save,Numberof_F/4,4,repeat_number3,plot_number);

train_set3 = 1:8;
test_set3 = 9:16;
S_save_mix = S_save;
second_pre_mix = second_term_test3;
dF2_pre_mix = dF2_pre3;
F_pre_mix = F_pre3;
Ffree_mix = F_save;
F_mix = F_save;
t00_mix = t00;
Range_mix = range_save0;

F_mix_test = reshape(F_mix(test_set3,:)',[1,length(test_set3)*discreteT_length]);
Fpre_mix_test = reshape(F_pre_mix(test_set3,:)',[1,length(test_set3)*discreteT_length]);

dF2train_mix = reshape(dF2_train3(train_set3,:)',[1,length(train_set3)*discreteT_length]);
dF2_mix = reshape(dF2dt2_save(train_set3,:)',[1,length(train_set3)*discreteT_length])/t_scale^2;


%% calculate err
%% train err
N_train_step = length(dF2_step);
error_RAEtrain_step = sum(abs(dF2_step- dF2train_step))/sum(abs(dF2_step));
error_RSEtrain_step = sum((dF2_step- dF2train_step).^2)/sum((dF2_step).^2);

N_train_trapz= length(dF2_mix);
error_RAEtrain_trapz= sum(abs(dF2_mix- dF2train_mix))/sum(abs(dF2_mix));
error_RSEtrain_trapz= sum((dF2_mix- dF2train_mix).^2)/sum((dF2_mix).^2);

N_train_mix = length(dF2_mix);
error_RAEtrain_mix = sum(abs(dF2_mix- dF2train_mix))/sum(abs(dF2_mix));
error_RSEtrain_mix = sum((dF2_mix- dF2train_mix).^2)/sum((dF2_mix).^2);



%% test err
N_test_step = length(F_step_test);
error_MAEtest_step = 1/N_test_step * sum(abs(F_step_test- Fpre_step_test)./abs(F_step_test));
error_MSEtest_step = 1/N_test_step * sum((F_step_test- Fpre_step_test).^2./(F_step_test).^2);
error_R2test_step = 1- sum((F_step_test- Fpre_step_test).^2)/sum((F_step_test-mean(F_step_test)).^2);

N_test_trapz = length(F_trapz_test);
error_MAEtest_trapz = 1/N_test_trapz * sum(abs(F_trapz_test- Fpre_trapz_test)./abs(F_trapz_test));
error_MSEtest_trapz = 1/N_test_trapz * sum((F_trapz_test- Fpre_trapz_test).^2./(F_trapz_test).^2);
error_R2test_trapz = 1- sum((F_trapz_test- Fpre_trapz_test).^2)/sum((F_trapz_test-mean(F_trapz_test)).^2);

N_test_mix = length(F_mix_test);
error_MAEtest_mix = 1/N_test_mix * sum(abs(F_mix_test- Fpre_mix_test)./abs(F_mix_test));
error_MSEtest_mix = 1/N_test_mix * sum((F_mix_test- Fpre_mix_test).^2./(F_mix_test).^2);
error_R2test_mix = 1- sum((F_mix_test- Fpre_mix_test).^2)/sum((F_mix_test-mean(F_mix_test)).^2);

%% err
err_all = [[N_train_step,error_RAEtrain_step,error_RSEtrain_step,0];
    [N_train_trapz,error_RAEtrain_trapz,error_RSEtrain_trapz,0];
    [N_train_mix,error_RAEtrain_mix,error_RSEtrain_mix,0];
    [N_test_step,error_MAEtest_step,error_MSEtest_step,error_R2test_step];
    [N_test_trapz,error_MAEtest_trapz,error_MSEtest_trapz,error_R2test_trapz]
    [N_test_mix,error_MAEtest_mix,error_MSEtest_mix,error_R2test_mix]]';

%% save
save(savePath,'dt00','t_scale','err_all','ep_number');
save(savePath,'test_set1','test_error1','err1','t00_step','Ffree_step','F_step','F_pre_step','S_save_step',...
    'Range_step','-append');
save(savePath,'test_set2','test_error2','err2','t00_trapz','Ffree_trapz','F_trapz','F_pre_trapz','S_save_trapz',...
    'Range_trapz','-append');
save(savePath,'test_set3','test_error3','err3','t00_mix','Ffree_mix','F_mix','F_pre_mix','S_save_mix',...
    'Range_mix','second_pre_mix','dF2_pre_mix','-append');


