function [err2,F_pre1,dF2_pre1,second_term_test1,dF2_train1] = fit2d_conti1(simu2,discreteT_length,Numberof_F,t00,t_scale,F_save,dFdt_save,dF2dt2_save,ds2dt2_save,Fi,Fj,repeat_number,plot_number)

dF2_pre1 = reshape(simu2(:,1)/t_scale^2,[discreteT_length Numberof_F])';
dF_pre1 = reshape(simu2(:,2)/t_scale,[discreteT_length Numberof_F])';
F_pre1 = reshape(simu2(:,3),[discreteT_length Numberof_F])';

% g3_test1 = reshape(simu2(:,4),[discreteT_length Numberof_F])';
g4_test1 = reshape(simu2(:,5),[discreteT_length Numberof_F])';
second_term_test1 = g4_test1.*ds2dt2_save/t_scale^2;

g3_train1 = reshape(simu2(:,6),[discreteT_length Numberof_F])';
g4_train1 = reshape(simu2(:,7),[discreteT_length Numberof_F])';
dF2_train1 = g3_train1+g4_train1.*ds2dt2_save;

second_term_train1 = g4_train1.*ds2dt2_save/t_scale^2;
dF2_train1 = dF2_train1/t_scale^2;

dF2dt2_save = dF2dt2_save/t_scale^2;
dFdt_save = dFdt_save/t_scale;

error_df2train = zeros(1,Numberof_F);
error_df2 = zeros(1,Numberof_F);
error_df = zeros(1,Numberof_F);
error_f = zeros(1,Numberof_F);

%% err
for i = 1:Numberof_F
    
    error_df2train(i) = sum((dF2dt2_save(i,:)- dF2_train1(i,:)).^2)/sum(dF2dt2_save(i,:).^2);
    error_df2(i) = sum((dF2dt2_save(i,:)- dF2_pre1(i,:)).^2)/sum(dF2dt2_save(i,:).^2);    
    error_df(i) = sum((dFdt_save(i,:)- dF_pre1(i,:)).^2)/sum(dFdt_save(i,:).^2);
    error_f(i) = mean((F_save(i,:)- F_pre1(i,:)).^2./F_save(i,:).^2);
end

err2 = [error_df2train;error_df2;error_df;error_f];

%% train
t1 = t00;
if plot_number == 1
    for i = 1:Numberof_F    
        figure(10*repeat_number+1);
        subplot(Fi,Fj,i);
        plot(t1,dF2dt2_save(i,:),'-','linewidth',1);hold on
        plot(t1,dF2_train1(i,:),'--','linewidth',1);hold on
        plot(t1,second_term_train1(i,:),':','linewidth',1);hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$d^2f/dt^2$','Interpreter','latex');
    end
    legend('Ref','Train 2d','$G_4s''$','Interpreter','latex');
    
    for i = 1:Numberof_F
        figure(10*repeat_number+2);
        subplot(Fi,Fj,i);
        plot(t1,dF2dt2_save(i,:),'-','linewidth',1);hold on
        plot(t1,dF2_pre1(i,:),'--','linewidth',1);hold on
        plot(t1,second_term_test1(i,:),':','linewidth',1);hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$d^2f/dt^2$','Interpreter','latex');
    end
    legend('Ref','Pre 2d','$G_4s''$','Interpreter','latex');
    
    % for i = 1:Numberof_F
    %     figure(10*repeat_number+3);
    %     subplot(Fi,Fj,i);
    %     plot(t1,dF2dt2_save(i,:),'-','linewidth',1);hold on
    %     plot(t1,dF2_pre1(i,:),'--','linewidth',1);hold on
    %     xlabel('$t$','Interpreter','latex');
    %     ylabel('$d^2f/dt^2$','Interpreter','latex');
    % end
    for i = 1:Numberof_F
        figure(10*repeat_number+4);
        subplot(Fi,Fj,i);
        plot(t1,F_save(i,:),'-','linewidth',1);hold on
        plot(t1,F_pre1(i,:),'--','linewidth',1);hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$f$','Interpreter','latex');  
    end
    legend('Ref','Pre 2d','Interpreter','latex');
end
end