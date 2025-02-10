function [err1,F_pre1,dF2_pre1,dF2_train1] = fit2d_disconti1(simu2,Length2,Numberof_F,t_scale,T_save2,F_save2,dFdt_save2,dF2dt2_save2,ds2dt2_save2,Fi,Fj,repeat_number,plot_number)

dF2_pre1 = simu2(:,1)'/t_scale^2;
dF_pre1 = simu2(:,2)'/t_scale;
F_pre1 = simu2(:,3)';

g3_pre1 = simu2(:,6)';
g4_pre1 = simu2(:,7)';

dF2_train1 = g3_pre1+g4_pre1.*ds2dt2_save2;
dF2_train1 = dF2_train1/t_scale^2;

dF2dt2_save2 = dF2dt2_save2/t_scale^2;
dFdt_save2 = dFdt_save2/t_scale;
T_save2 = T_save2*t_scale;

error_df2train = zeros(1,Numberof_F);
error_df2 = zeros(1,Numberof_F);
error_df = zeros(1,Numberof_F);
error_f = zeros(1,Numberof_F);

%%
for i = 1:Numberof_F
    num_t = Length2(i)+1:Length2(i+1);
    tX = T_save2(num_t);
    dF2dt2X = dF2dt2_save2(num_t);
    dFdtX = dFdt_save2(num_t);
    FX = F_save2(num_t);
    
    dF2_train1X = dF2_train1(num_t);
    dF2_pre1X = dF2_pre1(num_t);
    dF_pre1X = dF_pre1(num_t);
    F_pre1X = F_pre1(num_t);
    
    error_df2train(i) = sum((dF2dt2X- dF2_train1X).^2)/sum(dF2dt2X.^2);
    error_df2(i) = sum((dF2dt2X- dF2_pre1X).^2)/sum(dF2dt2X.^2);
    error_df(i) = sum((dFdtX- dF_pre1X).^2)/sum(dFdtX.^2);
    error_f(i) = mean((FX- F_pre1X).^2./FX.^2);
%%  
    if plot_number == 1
        figure(10*repeat_number+1);
        subplot(Fi,Fj,i);
        plot(tX,dF2dt2X,'b.');hold on
        plot(tX,dF2_train1X,'.');hold on
         xlabel('$t$','Interpreter','latex');
        ylabel('$d^2f/dt^2$','Interpreter','latex');
        figure(10*repeat_number+2);
        subplot(Fi,Fj,i);
        plot(tX,dF2dt2X,'.');hold on
        plot(tX,dF2_pre1X,'.');hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$d^2f/dt^2$','Interpreter','latex');
    %     figure(10*repeat_number+3);
    %     subplot(Fi,Fj,i);
    %     plot(tX,dFdtX,'.','linewidth',1);hold on
    %     plot(tX,dF_pre1X,'.');hold on
    %     xlabel('$t$','Interpreter','latex');
    %     ylabel('$df/dt$','Interpreter','latex'); 
        figure(10*repeat_number+4);
        subplot(Fi,Fj,i);
        plot(tX,FX,'.','linewidth',1);hold on
        plot(tX,F_pre1X,'.','linewidth',1);hold on
        xlabel('$t$','Interpreter','latex');
        ylabel('$f$','Interpreter','latex');
%         legend('Ref','Learned','Interpreter','latex');
    end
end
err1 = [error_df2train;error_df2;error_df;error_f];
end