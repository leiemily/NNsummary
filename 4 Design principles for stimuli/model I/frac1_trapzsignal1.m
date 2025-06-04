clear;
%% trapz signal
%% 条数，长度，Tp
Numberof_F = 16;
discreteT_length = 1023;
Numberof_Tp = 6;
% dt00 = 0.2;
load('frac1-1tp-mix1.mat','t_scale','dt00');

%% s setting
load('frac1-1tp-mix1.mat','range_save0');

step_maxvalue = 0.5;
step_minlen = 15*dt00;

slope_max = (0.1+max(abs(range_save0(Numberof_F+1,7:8))))/t_scale;
slope_min = 0.05/t_scale;

S_slope = zeros(Numberof_F,Numberof_Tp-1);
S_save = zeros(Numberof_F,discreteT_length);

for i = 1:Numberof_F
    [si,ss_slope,ss_steplen] = linear_s(dt00,discreteT_length,Numberof_Tp,step_maxvalue,step_minlen,slope_max,slope_min);
    S_slope(i,:) = ss_slope;
    S_save(i,:) = si;
end

%% figures
Fi = Numberof_F/2;
Fj = 2;
%% parameters
kR = 3.82e-2;
kBp = 3.25;

KR = 0.15;
Rt = 0.3;
Tt = 5/3;

m0 = 1;
Nr = 6;
alpha0 = 1.7;
KI = 18.2;
KA = 3000;

%% rescaling
t00 = 0:dt00:(discreteT_length-1)*dt00;
t10 = t00/t_scale;
dt10 = dt00/t_scale;

%% cut
Xk_cut1 = 2;
Xk_cut2 = 1;
%% save space
F_save = zeros(Numberof_F,discreteT_length);
m_save = zeros(Numberof_F,discreteT_length);
Yp_save = zeros(Numberof_F,discreteT_length);


%% iteration
for i = 1:Numberof_F
    %% linear signal s(t)
    s = S_save(i,:);

   %% save space
    m = zeros(1,discreteT_length);
    A = zeros(1,discreteT_length);
    E = zeros(1,discreteT_length);

    R = zeros(1,discreteT_length);
    Bp = zeros(1,discreteT_length);
    Yp = zeros(1,discreteT_length);

    %% initial space
%     m_initial = 0.71;
    m_initial = 0.7047;
    m(1) = m_initial;
    E(1) = - alpha0 *(m(1) - m0) + log((1 + s(1)/KI)/(1 + s(1)/KA));
    A(1) = (1 + exp(Nr * E(1)))^(-1);
    R(1) = Rt / (1 + KR * Tt * (1 - A(1)));

    x0 = [1,1,1,1,1,1];
    %% iteration t
    for j = 1:discreteT_length-1
        c = A(j);
        fun = @(x)paramfun1(x,c);
        S = fsolve(fun,x0);
        x0 = S;

        Bp(j) = S(5);
        Yp(j) = S(6);

        m(j+1)  = m(j) + dt00 * (kR * R(j) * (1 - A(j)) - kBp * Bp(j) * A(j));
        E(j+1) = - alpha0 *(m(j+1) - m0) + log((1 + s(j+1)/KI)/(1 + s(j+1)/KA));
        A(j+1) = (1 + exp(Nr * E(j+1)))^(-1);
        R(j+1) = Rt / (1 + KR * Tt * (1 - A(j+1)));

    end
    c = A(end);
    fun = @(x)paramfun1(x,c);
    S = fsolve(fun,x0);
    Yp(end) = S(6);

    F = cellmovemodel_3(Yp);
    %% update space
    F_save(i,:) = F;
    m_save(i,:) = m;
    Yp_save(i,:) = Yp;
end
%% 2系列代表筛选出的数据集
T_save2 = [];
S_save2 = [];
F_save2 = [];
dFdt_save2 = [];
dF2dt2_save2 = [];
dsdt_save2 = [];
ds2dt2_save2 = [];

Length2 = zeros(1,Numberof_F+1);
Length2lls = zeros(Numberof_F,2*Numberof_Tp);
range_save0 = zeros(Numberof_F+1,10);
%%
for i = 1:Numberof_F
    m = m_save(i,:);
    Yp = Yp_save(i,:);
    F = F_save(i,:);
    s = S_save(i,:);
    %% gradients 修改了求导方式！！！！！！！！！！！！
    dsdt = gradient1_s(s,dt10);
    ds2dt2 = gradient2_s(s,dt10);

    lls = [0,find(abs(ds2dt2)>1e-2),discreteT_length];
    %% filtering and smooth
    X = [];
    dFdt = [];
    dF2dt2 = [];

    for k = 1:length(lls)-1
        Xk =  Xk_cut1 + lls(k):lls(k+1)-Xk_cut2;
        X = [X,Xk];
        Length2lls(i,k+1) = Length2lls(i,k) + length(Xk);
        dFdt = [dFdt, gradient1_f(F(Xk),dt10)];
        dF2dt2 = [dF2dt2, gradient2_f(F(Xk),dt10)];
    end
    %% update space    
    T_save2 = [T_save2,t10(X)];
    S_save2 = [S_save2,s(X)];
    F_save2 = [F_save2,F(X)];
    dFdt_save2 = [dFdt_save2,dFdt];
    dF2dt2_save2 = [dF2dt2_save2,dF2dt2];
    dsdt_save2 = [dsdt_save2, dsdt(X)];
    ds2dt2_save2 = [ds2dt2_save2, zeros(1,length(X))];

    Length2(i+1) = length(T_save2);
    %% save range
    range_save0(i,:) = [min(F),max(F),min(s(X)),max(s(X)),...
        min(dFdt),max(dFdt),min(dsdt(X)),max(dsdt(X)),...
        min(dF2dt2),max(dF2dt2)];
    %% figures
    figure(1);   
    subplot(Fi,Fj,i);
    plot(t00,m,'-','linewidth',2); hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$m$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,Yp_save(i,:),'--');hold on
    ylabel('$Yp$','Interpreter','latex','Fontsize',18);   
    

    figure(2);   
    subplot(Fi,Fj,i);
%     plot(t00,F_noise,'.'); hold on
    plot(t00,F,'-','linewidth',2); hold on
%     plot(t00(X),F_smooth,'.');hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$f$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,s,'--','linewidth',2);hold on
    ylabel('$s$','Interpreter','latex','Fontsize',18);

    figure(3);   
    subplot(Fi,Fj,i);
    plot(t10(X),dFdt,'.');hold on
%     plot(t10(X),dFdt_smooth,'.');hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$df/dt$','Interpreter','latex','Fontsize',18);

    figure(4);   
    subplot(Fi,Fj,i);
    plot(t10(X),dF2dt2,'.');hold on
%     plot(t10(X),dF2dt2_smooth,'.');
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$d^2f/dt^2$','Interpreter','latex','Fontsize',18); 
end
%%
range_save0(Numberof_F+1,:) = [min(range_save0(:,1)),max(range_save0(:,2)),min(range_save0(:,3)),max(range_save0(:,4)),...
    min(range_save0(:,5)),max(range_save0(:,6)),min(range_save0(:,7)),max(range_save0(:,8)),...
    min(range_save0(:,9)),max(range_save0(:,10))];
%%
save('frac1-1tp-trapz1.mat','t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'm_save','Yp_save','Length2','Length2lls','T_save2',...
    'S_save','S_save2','dsdt_save2','ds2dt2_save2',...
    'F_save','F_save2','dFdt_save2','dF2dt2_save2',...
    'range_save0');