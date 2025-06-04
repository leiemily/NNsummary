clear;
%% main
j_set = 10:2:24;
plot_number = 0;%plot when ==1 else do not plot
test_set3 = 1:4;

err = zeros(4*length(j_set),5);
%% paths
folderPath_varieddata = '.\dataset1\';
folderPath_NNdata = '.\dataset2\';

%%
for kj = 1:length(j_set)
    j = j_set(kj);
    trainsetName0 = ['frac2-',num2str(j),'tps-save1d-free'];
    load(fullfile(folderPath_NNdata, trainsetName0));

    testsetName0 = ['frac2-',num2str(j),'tps-test1d-free'];
    savePath = fullfile(folderPath_NNdata, testsetName0);
    if kj == 1
        datasetName3 = ['frac2-',num2str(j),'tpb-expmix1'];
    else
        datasetName3 = ['frac2-',num2str(j),'tps-expmix1'];
    end
    
    load(fullfile(folderPath_varieddata, datasetName3));
    simu_expcos = Y_all3;
    repeat_number3 = 3;

    [err3,F_pre3,dF_pre3,second_term_test3,dF_train3] = fit1d_conti1(simu_expcos,discreteT_length,Numberof_F,t00,t_scale,...
        F_save,dFdt_save,dsdt_save,Numberof_F/2,2,repeat_number3,plot_number);

    
    S_save_expcos = S_save;
    second_pre_expcos = second_term_test3;
    dF_pre_expcos = dF_pre3;
    F_pre_expcos = F_pre3;
    F_expcos = F_save;
    t00_expcos = t00;

    F_expcos_test = reshape(F_expcos(test_set3,:)',[1,length(test_set3)*discreteT_length]);
    Fpre_expcos_test = reshape(F_pre_expcos(test_set3,:)',[1,length(test_set3)*discreteT_length]);
    Range_expcos = range_save0;

    %% test err
    N_test_expcos = length(F_expcos_test);
    error_MAEtest_expcos = 1/N_test_expcos * sum(abs(F_expcos_test- Fpre_expcos_test)./abs(F_expcos_test));
    error_MSEtest_expcos = 1/N_test_expcos * sum((F_expcos_test- Fpre_expcos_test).^2./(F_expcos_test).^2);
    error_R2test_expcos = 1- sum((F_expcos_test- Fpre_expcos_test).^2)/sum((F_expcos_test-mean(F_expcos_test)).^2);

    err_all = [N_test_expcos,error_MAEtest_expcos,error_MSEtest_expcos,error_R2test_expcos]';
    err(4*(kj-1)+1:4*kj,:) = [[zeros(1,4);err3],err_all];
    
    %% save
    save(savePath,'dt00','t_scale','err_all','ep_number');
    save(savePath,'test_set3','test_error3','err3','t00_expcos','F_expcos','F_pre_expcos','S_save_expcos',...
        'Range_expcos','second_pre_expcos','dF_pre_expcos','-append');
end
