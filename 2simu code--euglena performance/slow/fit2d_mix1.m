%%
simu2 = Y_all2;

dF2_pre1 = reshape(simu2(:,1)/t_scale^2,[discreteT_length Numberof_F])';
dF_pre1 = reshape(simu2(:,2)/t_scale,[discreteT_length Numberof_F])';
F_pre1 = reshape(simu2(:,3),[discreteT_length Numberof_F])';

g3_pre1 = reshape(simu2(:,6),[discreteT_length Numberof_F])';
g4_pre1 = reshape(simu2(:,7),[discreteT_length Numberof_F])';
dF2_train1 = g3_pre1+g4_pre1.*ds2dt2_save;
second_term = g4_pre1.*ds2dt2_save/t_scale^2;
dF2_train1 = dF2_train1/t_scale^2;

dF2dt2smooth_save = dF2dt2smooth_save/t_scale^2;
dFdtsmooth_save = dFdtsmooth_save/t_scale;
ds2dt2_save = ds2dt2_save/t_scale^2;
error_df2train = zeros(1,Numberof_F);
error_df2 = zeros(1,Numberof_F);
error_df = zeros(1,Numberof_F);
error_f = zeros(1,Numberof_F);

%% err
for i = 1:Numberof_F
    
    error_df2train(i) = sum((dF2dt2smooth_save(i,:)- dF2_train1(i,:)).^2)/sum(dF2dt2smooth_save(i,:).^2);
    error_df2(i) = sum((dF2dt2smooth_save(i,:)- dF2_pre1(i,:)).^2)/sum(dF2dt2smooth_save(i,:).^2);    
    error_df(i) = sum((dFdtsmooth_save(i,:)- dF_pre1(i,:)).^2)/sum(dFdtsmooth_save(i,:).^2);
%     error_f(i) = sum((Fsmooth_save(i,:)- F_pre1(i,:)).^2)/sum(Fsmooth_save(i,:).^2);
    error_f(i) = mean((Fsmooth_save(i,:)- F_pre1(i,:)).^2./Fsmooth_save(i,:).^2);
end

err2 = [error_df2train;error_df2;error_df;error_f];

%% train
Fi = 3;
Fj = 3;
t1 = t00;
for i = 1:Numberof_F
    figure(10);
    subplot(Fi,Fj,i);
% %     plot(t1,dF2dt2smooth_save(i,:),'k-','linewidth',1);hold on
%     plot(t1,g3_pre1(i,:),'-','linewidth',1);hold on
%     ylabel('$G_3$','Interpreter','latex');
%     yyaxis right
%     plot(t1,Fs_save(i,:),'-','linewidth',1);hold on
    plot(t1,second_term(i,:),'-','linewidth',1);hold on
    plot(t1,dF2_train1(i,:),'--','linewidth',1);hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$G_4s''$','Interpreter','latex');
    figure(11);
    subplot(Fi,Fj,i);
    plot(t1,dF2dt2smooth_save(i,:),'-','linewidth',1);hold on
    plot(t1,dF2_train1(i,:),'--','linewidth',1);hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$d^2f/dt^2$','Interpreter','latex');
end

for i = 1:Numberof_F
    figure(12);
    subplot(Fi,Fj,i);
    plot(t1,dF2dt2smooth_save(i,:),'-','linewidth',1);hold on
    plot(t1,dF2_pre1(i,:),'--','linewidth',1);hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$d^2f/dt^2$','Interpreter','latex');
end
for i = 1:Numberof_F
    figure(14);
    subplot(Fi,Fj,i);
    plot(t1,Fsmooth_save(i,:),'-','linewidth',1);hold on
    plot(t1,F_pre1(i,:),'--','linewidth',1);hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$f$','Interpreter','latex');  
end
legend('Ref','Pre 2d','Interpreter','latex');