%%
simu2 = Y_all2;

dF_pre1 = reshape(simu2(:,1),[discreteT_length Numberof_F])'/t_scale;
F_pre1 = reshape(simu2(:,2),[discreteT_length Numberof_F])';
g1_pre1 = reshape(simu2(:,5),[discreteT_length Numberof_F])';
g2_pre1 = reshape(simu2(:,6),[discreteT_length Numberof_F])';
dF_train1 = g1_pre1 + g2_pre1.*dsdt_save;
dF_train1 = dF_train1/t_scale;

dFdtsmooth_save = dFdtsmooth_save/t_scale;
dsdt_save = dsdt_save/t_scale;
error_dftrain = zeros(1,Numberof_F);
error_df = zeros(1,Numberof_F);
error_f = zeros(1,Numberof_F);
%% err
for i = 1:Numberof_F
   
    error_dftrain(i) = sum((dFdtsmooth_save(i,:)- dF_train1(i,:)).^2)/sum(dFdtsmooth_save(i,:).^2);
    error_df(i) = sum((dFdtsmooth_save(i,:)- dF_pre1(i,:)).^2)/sum(dFdtsmooth_save(i,:).^2);
%     error_f(i) = sum((Fsmooth_save(i,:)- F_pre1(i,:)).^2)/sum(Fsmooth_save(i,:).^2);
    error_f(i) = mean((Fsmooth_save(i,:)- F_pre1(i,:)).^2./Fsmooth_save(i,:).^2);
end

err2 = [error_dftrain;error_df;error_f];

%% train
Fi = 3;
Fj = 6;
t1 = t00;
for i = 1:Numberof_F
    figure(10);
    subplot(Fi,Fj,i);
%     plot(t1,Fs_save(i,:),'-','linewidth',1);hold on
    plot(t1,g2_pre1(i,:),'--','linewidth',1);hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$G_2$','Interpreter','latex');
    
    figure(11);
    subplot(Fi,Fj,i);
    plot(t1,dFdtsmooth_save(i,:),'-','linewidth',1);hold on
    plot(t1,dF_train1(i,:),'--','linewidth',1);hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$df/dt$','Interpreter','latex');
end
%% test
for i = 1:Numberof_F
    figure(13);
    subplot(Fi,Fj,i);
    plot(t1,dFdtsmooth_save(i,:),'-','linewidth',1);hold on
    plot(t1,dF_pre1(i,:),'--','linewidth',1);hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$df/dt$','Interpreter','latex');
    figure(14);
    subplot(Fi,Fj,i);
    plot(t1,Fsmooth_save(i,:),'-','linewidth',1);hold on
    plot(t1,F_pre1(i,:),'--','linewidth',1);hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$f$','Interpreter','latex');  
end
legend('Ref','Learned','Interpreter','latex');