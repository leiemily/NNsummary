clear;
%% main
j = 10;
plot_number = 1;%plot when ==1 else do not plot
%% paths
folderPath_varieddata = '.\dataset1\';
folderPath_NNdata = '.\dataset2\';

trainsetName0 = ['bias2-',num2str(j),'tp-save2d-noise'];
load(fullfile(folderPath_NNdata, trainsetName0));

testsetName0 = ['bias2-',num2str(j),'tp-test2d-noise'];
savePath = fullfile(folderPath_NNdata, testsetName0);

%% test errs
Fi = 2;
Fj = 2;
figure(1);
subplot(Fi,Fj,1);
plot(ep_number, test_error1,'o-');
set(gca, 'YScale', 'log');

subplot(Fi,Fj,2);
plot(ep_number, test_error2,'o-');
set(gca, 'YScale', 'log');

subplot(Fi,Fj,3);
plot(ep_number, test_error3,'o-');
set(gca, 'YScale', 'log');

subplot(Fi,Fj,4);
plot(ep_number, test_error4,'o-');
set(gca, 'YScale', 'log');

%%
datasetName1 = ['bias2-',num2str(j),'tp-step1'];
load(fullfile(folderPath_varieddata, datasetName1));
simu_step = Y_all1;
repeat_number1 = 1;
[err1,F_pre1,dF2_train1] = fit2d_disconti1(simu_step,Length2,Numberof_F,t_scale,T_save2,...
    Fsmooth_save2,dFdtsmooth_save2,dF2dt2smooth_save2,ds2dt2_save2,Numberof_F/2,2,repeat_number1,plot_number);

test_set1 = 7:10;
S_save_step = S_save;
F_pre_step = F_pre1;
Fnoise_step = Fnoise_save;
F_step = Fsmooth_save2;
t00_step = t00;
Range_step = range_save;

num_test_set1 = Length2(test_set1(1))+1:Length2(test_set1(end)+1);
F_step_test = F_step(num_test_set1);
Fpre_step_test = F_pre_step(num_test_set1);

train_set1 = 1:6;
num_train_set1 = Length2(train_set1(1))+1:Length2(train_set1(end)+1);
dF2train_step = dF2_train1(num_train_set1);
dF2_step = dF2dt2smooth_save2(num_train_set1)/t_scale^2;

%%
datasetName2 = ['bias2-',num2str(j),'tp-mix1'];
load(fullfile(folderPath_varieddata, datasetName2));
simu_mix = Y_all2;
repeat_number2 = 2;
[err2,F_pre1,dF2_train1] = fit2d_conti1(simu_mix,discreteT_length,Numberof_F,t00,t_scale,...
    Fsmooth_save,dFdtsmooth_save,dF2dt2smooth_save,ds2dt2_save,Numberof_F/2,2,repeat_number2,plot_number);


test_set2 = 9:10;
S_save_mix = S_save;
F_pre_mix = F_pre1;
Fnoise_mix = Fnoise_save;
F_mix = Fsmooth_save;
t00_mix = t00;
Range_mix = range_save;

F_mix_test = reshape(F_mix(test_set2,:)',[1,length(test_set2)*discreteT_length]);
Fpre_mix_test = reshape(F_pre_mix(test_set2,:)',[1,length(test_set2)*discreteT_length]);

train_set2 = 1:8;
dF2train_mix = reshape(dF2_train1(train_set2,:)',[1,length(train_set2)*discreteT_length]);
dF2_mix = reshape(dF2dt2smooth_save(train_set2,:)',[1,length(train_set2)*discreteT_length])/t_scale^2;

%%
datasetName3 = ['bias2-',num2str(j),'tp-expmix1'];
load(fullfile(folderPath_varieddata, datasetName3));
simu_expcos = Y_all3;
repeat_number3 = 3;
[err3,F_pre1,dF2_train1] = fit2d_conti1(simu_expcos,discreteT_length,Numberof_F,t00,t_scale,...
    Fsmooth_save,dFdtsmooth_save,dF2dt2smooth_save,ds2dt2_save,Numberof_F/2,2,repeat_number3,plot_number);

test_set3 = 1:4;
S_save_expcos = S_save;
F_pre_expcos = F_pre1;
Fnoise_expcos = Fnoise_save;
F_expcos = Fsmooth_save;
t00_expcos = t00;

F_expcos_test = reshape(F_expcos(test_set3,:)',[1,length(test_set3)*discreteT_length]);
Fpre_expcos_test = reshape(F_pre_expcos(test_set3,:)',[1,length(test_set3)*discreteT_length]);
Range_expcos = range_save;

%%
datasetName4 = ['bias2-',num2str(j),'tp-gauss1'];
load(fullfile(folderPath_varieddata, datasetName4));
simu_gauss = Y_all4;
repeat_number4 = 4;
[err4,F_pre1,dF2_train1] = fit2d_conti1(simu_gauss,discreteT_length,Numberof_F,t00,t_scale,...
    Fsmooth_save,dFdtsmooth_save,dF2dt2smooth_save,ds2dt2_save,Numberof_F/2,2,repeat_number4,plot_number);

test_set4 = 1:4;
S_save_gauss = S_save;
F_pre_gauss = F_pre1;
Fnoise_gauss = Fnoise_save;
F_gauss = Fsmooth_save;
t00_gauss = t00;

F_gauss_test = reshape(F_gauss(test_set4,:)',[1,length(test_set4)*discreteT_length]);
Fpre_gauss_test = reshape(F_pre_gauss(test_set4,:)',[1,length(test_set4)*discreteT_length]);
Range_gauss = range_save;


%% calculate err
%% train err
N_train_step = length(dF2_step);
N_train_mix = length(dF2_mix);
train_all = [dF2_step,dF2_mix];
trainpre_all = [dF2train_step,dF2train_mix];
N_train_all = length(train_all);


error_RAEtrain_step = sum(abs(dF2_step- dF2train_step))/sum(abs(dF2_step));
error_RSEtrain_step = sum((dF2_step- dF2train_step).^2)/sum((dF2_step).^2);

error_RAEtrain_mix = sum(abs(dF2_mix- dF2train_mix))/sum(abs(dF2_mix));
error_RSEtrain_mix = sum((dF2_mix- dF2train_mix).^2)/sum((dF2_mix).^2);

error_RAEtrain_all = sum(abs(train_all- trainpre_all))/sum(abs(train_all));
error_RSEtrain_all = sum((train_all- trainpre_all).^2)/sum((train_all).^2);

%% test err
N_test_step = length(F_step_test);
error_MAEtest_step = 1/N_test_step * sum(abs(F_step_test- Fpre_step_test)./abs(F_step_test));
error_MSEtest_step = 1/N_test_step * sum((F_step_test- Fpre_step_test).^2./(F_step_test).^2);
error_R2test_step = 1- sum((F_step_test- Fpre_step_test).^2)/sum((F_step_test-mean(F_step_test)).^2);

N_test_mix = length(F_mix_test);
error_MAEtest_mix = 1/N_test_mix * sum(abs(F_mix_test- Fpre_mix_test)./abs(F_mix_test));
error_MSEtest_mix = 1/N_test_mix * sum((F_mix_test- Fpre_mix_test).^2./(F_mix_test).^2);
error_R2test_mix = 1- sum((F_mix_test- Fpre_mix_test).^2)/sum((F_mix_test-mean(F_mix_test)).^2);

N_test_expcos = length(F_expcos_test);
error_MAEtest_expcos = 1/N_test_expcos * sum(abs(F_expcos_test- Fpre_expcos_test)./abs(F_expcos_test));
error_MSEtest_expcos = 1/N_test_expcos * sum((F_expcos_test- Fpre_expcos_test).^2./(F_expcos_test).^2);
error_R2test_expcos = 1- sum((F_expcos_test- Fpre_expcos_test).^2)/sum((F_expcos_test-mean(F_expcos_test)).^2);

N_test_gauss = length(F_gauss_test);
error_MAEtest_gauss = 1/N_test_gauss * sum(abs(F_gauss_test- Fpre_gauss_test)./abs(F_gauss_test));
error_MSEtest_gauss = 1/N_test_gauss * sum((F_gauss_test- Fpre_gauss_test).^2./(F_gauss_test).^2);
error_R2test_gauss = 1- sum((F_gauss_test- Fpre_gauss_test).^2)/sum((F_gauss_test -mean(F_gauss_test)).^2);

%% err
err_all = [[N_train_step,error_RAEtrain_step,error_RSEtrain_step,0];
    [N_train_mix,error_RAEtrain_mix,error_RSEtrain_mix,0];
    [N_train_all,error_RAEtrain_all,error_RSEtrain_all,0];
    [N_test_mix,error_MAEtest_mix,error_MSEtest_mix,error_R2test_mix];
    [N_test_expcos,error_MAEtest_expcos,error_MSEtest_expcos,error_R2test_expcos];
    [N_test_gauss,error_MAEtest_gauss,error_MSEtest_gauss,error_R2test_gauss];
    [N_test_step,error_MAEtest_step,error_MSEtest_step,error_R2test_step]]';

%% save
save(savePath,'dt00','t_scale','err_all','ep_number');
save(savePath,'test_set1','test_error1','err1','t00_step','Fnoise_step','F_step','F_pre_step','S_save_step','Range_step','-append');
save(savePath,'test_set2','test_error2','err2','t00_mix','Fnoise_mix','F_mix','F_pre_mix','S_save_mix','Range_mix','-append');
save(savePath,'test_set3','test_error3','err3','t00_expcos','Fnoise_expcos','F_expcos','F_pre_expcos','S_save_expcos','Range_expcos','-append');
save(savePath,'test_set4','test_error4','err4','t00_gauss','Fnoise_gauss','F_gauss','F_pre_gauss','S_save_gauss','Range_gauss','-append');
