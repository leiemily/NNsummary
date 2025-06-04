function [err1,F_pre1,dF_pre1,dF_train1] = fit1d_disconti1(simu2,Length2,Numberof_F,t_scale,T_save2,F_save2,dFdt_save2,dsdt_save2,Fi,Fj,repeat_number,plot_number)

dF_pre1 = simu2(:,1)'/t_scale;
F_pre1 = simu2(:,2)';

g1_pre1 = simu2(:,5)';
g2_pre1 = simu2(:,6)';

dF_train1 = g1_pre1 + g2_pre1.*dsdt_save2;
dF_train1 = dF_train1/t_scale;

dFdt_save2 = dFdt_save2/t_scale;
T_save2 = T_save2*t_scale;

error_dftrain = zeros(1,Numberof_F);
error_df = zeros(1,Numberof_F);
error_f = zeros(1,Numberof_F);

%%
for i = 1:Numberof_F
    num_t = Length2(i)+1:Length2(i+1);
    tX = T_save2(num_t);
    dFdtX = dFdt_save2(num_t);
    FX = F_save2(num_t);
    
    dF_train1X = dF_train1(num_t);
    dF_pre1X = dF_pre1(num_t);
    F_pre1X = F_pre1(num_t);
    error_dftrain(i) = sum((dFdtX- dF_train1X).^2)/sum(dFdtX.^2);
    error_df(i) = sum((dFdtX- dF_pre1X).^2)/sum(dFdtX.^2);
    error_f(i) = mean((FX- F_pre1X).^2./FX.^2);
    %%
    if plot_number == 1
        figure(10*repeat_number+1);
        subplot(Fi,Fj,i);
        plot(tX,dFdtX,'.');hold on
        plot(tX,dF_train1X,'.');hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$df/dt$','Interpreter','latex');

%         figure(10*repeat_number+2);
%         subplot(Fi,Fj,i);
%         plot(tX,dFdtX,'.','linewidth',1);hold on
%         plot(tX,dF_pre1X,'.');hold on
%         xlabel('$t$','Interpreter','latex');
%         ylabel('$df/dt$','Interpreter','latex');

        figure(10*repeat_number+3);
        subplot(Fi,Fj,i);
        plot(tX,FX,'.','linewidth',1);hold on
        plot(tX,F_pre1X,'.','linewidth',1);hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$f$','Interpreter','latex'); 
        legend('Ref','Learned','Interpreter','latex');
    end
end
err1 = [error_dftrain;error_df;error_f];
end