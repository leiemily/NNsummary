function [err2,F_pre1,dF_train1] = fit1d_conti1(simu2,discreteT_length,Numberof_F,t00,t_scale,F_save,dFdt_save,dsdt_save,Fi,Fj,repeat_number,plot_number)

dF_pre1 = reshape(simu2(:,1),[discreteT_length Numberof_F])'/t_scale;
F_pre1 = reshape(simu2(:,2),[discreteT_length Numberof_F])';

g1_pre1 = reshape(simu2(:,5),[discreteT_length Numberof_F])';
g2_pre1 = reshape(simu2(:,6),[discreteT_length Numberof_F])';

dF_train1 = g1_pre1 + g2_pre1.*dsdt_save;
dF_train1 = dF_train1/t_scale;

dFdt_save = dFdt_save/t_scale;
error_dftrain = zeros(1,Numberof_F);
error_df = zeros(1,Numberof_F);
error_f = zeros(1,Numberof_F);
%% err
for i = 1:Numberof_F  
    error_dftrain(i) = sum((dFdt_save(i,:)- dF_train1(i,:)).^2)/sum(dFdt_save(i,:).^2);
    error_df(i) = sum((dFdt_save(i,:)- dF_pre1(i,:)).^2)/sum(dFdt_save(i,:).^2);
    error_f(i) = mean((F_save(i,:)- F_pre1(i,:)).^2./F_save(i,:).^2);
end

err2 = [error_dftrain;error_df;error_f];

%% train
t1 = t00;
if plot_number == 1
    for i = 1:Numberof_F
%         figure(10*repeat_number+1);
%         subplot(Fi,Fj,i);
%     %     plot(t1,Fs_save(i,:),'-','linewidth',1);hold on
%         plot(t1,g2_pre1(i,:),'--','linewidth',1);hold on
%         xlabel('$t$','Interpreter','latex');
%         ylabel('$G_2$','Interpreter','latex');

        figure(10*repeat_number+2);
        subplot(Fi,Fj,i);
        plot(t1,dFdt_save(i,:),'-','linewidth',1);hold on
        plot(t1,dF_train1(i,:),'--','linewidth',1);hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$df/dt$','Interpreter','latex');
    end
    %% test
    for i = 1:Numberof_F
%         figure(10*repeat_number+3);
%         subplot(Fi,Fj,i);
%         plot(t1,dFdt_save(i,:),'-','linewidth',1);hold on
%         plot(t1,dF_pre1(i,:),'--','linewidth',1);hold on
%         xlabel('$t$','Interpreter','latex');
%         ylabel('$df/dt$','Interpreter','latex');
        figure(10*repeat_number+4);
        subplot(Fi,Fj,i);
        plot(t1,F_save(i,:),'-','linewidth',1);hold on
        plot(t1,F_pre1(i,:),'--','linewidth',1);hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$f$','Interpreter','latex');  
    end
    legend('Ref','Learned','Interpreter','latex');
end
end