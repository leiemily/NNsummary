clear;
%% step signals
%% 条数，长度，Tp
Numberof_F = 10;
discreteT_length = 1002;
Numberof_Tp = 10;
load('bias2-10tpn-mix1.mat','t_scale','dt00');

%% noise level
% nl_snr = unifrnd(20,22,1,Numberof_F);
nl_snr = 22*ones(1,Numberof_F);
snr_linear = 10.^(nl_snr/10);% 将SNR转换为线性比例

%% s setting
step_maxvalue = 2;
step_minlen = (discreteT_length-1)*dt00/15;

rng("shuffle");
S_para = unifrnd(0,step_maxvalue,Numberof_F, Numberof_Tp-1);
S_para = [zeros(Numberof_F,1),S_para];

Tp = zeros(Numberof_F,Numberof_Tp+1);
for i = 1:Numberof_F
    Tpi_diff = zeros(1,Numberof_Tp);
    while min(Tpi_diff) < step_minlen
        Tpi = randi([1,fix((discreteT_length-1)*dt00)],1, Numberof_Tp-1);
        Tpi = sort(Tpi,2);
        Tpi = [0,Tpi,(discreteT_length-1)*dt00];
        Tpi_diff = diff(Tpi);
    end
    Tp(i,:) = Tpi;
end

% load('bias2-10tp-step1.mat','S_para','Tp');

%% figures
Fi = Numberof_F/2;
Fj = 2;

%% parameters
m0 = 1;
Nr = 6;
alpha0 = 1.7;
KI = 18.2;
KA = 3000;

V_R = 0.01;
V_B = 0.02;
V_on = 0.02;
V_off = 0.0194;
FliM_max = 45;

K_A = 0.35;
P = 1e7;
C = 4.1;

aa = 30;
bb = 8;

%% rescaling
t00 = 0:dt00:(discreteT_length-1)*dt00;
t10 = t00/t_scale;
dt10 = dt00/t_scale;
%% cut
Xk_cut1 = 2;
Xk_cut2 = 1;
%% save space
S_save = zeros(Numberof_F,discreteT_length);%储存加噪s
F_save = zeros(Numberof_F,discreteT_length);%储存加噪s下的f
S0_save = zeros(Numberof_F,discreteT_length);%储存无噪s
F0_save = zeros(Numberof_F,discreteT_length);%储存无噪s下的f
m_save = zeros(Numberof_F,discreteT_length);
FliM_save = zeros(Numberof_F,discreteT_length);

%% iteration
%% noiseless
for i = 1:Numberof_F
    %% step signal s(t)
    ai = S_para(i,:);
    Tpi = Tp(i,:);
    s = zeros(1,discreteT_length);
    for k = 1:length(Tpi)-1
        s = s + ai(k)*(t00>=Tpi(k) & t00<Tpi(k+1));
    end
    s(end) = ai(end);
 
   %% save space
    m = zeros(1,discreteT_length);
    FliM = zeros(1,discreteT_length);
    A = zeros(1,discreteT_length);
    E = zeros(1,discreteT_length);

    %% initial space
    x0 = [1,34];
    c = s(1);
    fun = @(x)paramfun_m(x,c);
    x0 = fsolve(fun,x0);
    m_initial = x0(1);
    FliM_initial = x0(2);
    
    m(1) = m_initial;
    FliM(1) = FliM_initial;
    E(1) = - alpha0 *(m(1) - m0) + log((1 + s(1)/KI)/(1 + s(1)/KA));
    A(1) = (1 + exp (Nr * E(1))).^(-1);

    %% iteration
    for j = 1:discreteT_length-1
        m(j+1)  = m(j) + dt00 * (V_R*(1 - A(j))- V_B*A(j));
        FliM(j+1) = FliM(j) + dt00 * (V_on * (FliM_max - FliM(j)) - V_off * A(j)* FliM(j));

        E(j+1) = - alpha0 *(m(j+1) - m0) + log((1 + s(j+1)/KI)/(1 + s(j+1)/KA));
        A(j+1) = (1 + exp (Nr * E(j+1))).^(-1);
    end

    %% calculate F
    F = (1 + A/K_A).^((FliM-aa)*bb) ./ ((1 + A/K_A).^((FliM-aa)*bb) + P * (1 + A/(K_A*C)).^((FliM-aa)*bb));
    %% update space
    S0_save(i,:) = s;
    F0_save(i,:) = F;
end

%% noisy
p_signal = zeros(1,Numberof_F);
p_noise = zeros(1,Numberof_F);

for i = 1:Numberof_F
    %% step signal s(t)
    s_origin = S0_save(i,:);
    
    %% add  noise
    p_signal(i) = mean(s_origin.^2);
    p_noise(i) = p_signal(i) / snr_linear(i);% 计算噪声功率
    sigma = sqrt(p_noise(i));
    
    rng('default');
    seed_noise = rng;
    rng(seed_noise);
    noise_s = normrnd(0,sigma,size(s_origin));
    s = s_origin + noise_s;
    s_overnn = find(s<0);
    s(s_overnn) = s_origin(s_overnn);

   %% save space
    m = zeros(1,discreteT_length);
    FliM = zeros(1,discreteT_length);
    A = zeros(1,discreteT_length);
    E = zeros(1,discreteT_length);

    %% initial space
    x0 = [1,34];
    c = s(1);
    fun = @(x)paramfun_m(x,c);
    x0 = fsolve(fun,x0);
    m_initial = x0(1);
    FliM_initial = x0(2);
    
    m(1) = m_initial;
    FliM(1) = FliM_initial;
    E(1) = - alpha0 *(m(1) - m0) + log((1 + s(1)/KI)/(1 + s(1)/KA));
    A(1) = (1 + exp (Nr * E(1))).^(-1);

    %% iteration
    for j = 1:discreteT_length-1
        m(j+1)  = m(j) + dt00 * (V_R*(1 - A(j))- V_B*A(j));
        FliM(j+1) = FliM(j) + dt00 * (V_on * (FliM_max - FliM(j)) - V_off * A(j)* FliM(j));

        E(j+1) = - alpha0 *(m(j+1) - m0) + log((1 + s(j+1)/KI)/(1 + s(j+1)/KA));
        A(j+1) = (1 + exp (Nr * E(j+1))).^(-1);
    end

    %% calculate F
    F = (1 + A/K_A).^((FliM-aa)*bb) ./ ((1 + A/K_A).^((FliM-aa)*bb) + P * (1 + A/(K_A*C)).^((FliM-aa)*bb));
    %% update space
    S_save(i,:) = s;
    F_save(i,:) = F;
    m_save(i,:) = m;
    FliM_save(i,:) = FliM;
end

%% noise and smooth  noise and smooth  noise and smooth  
p_smooth = 0.005*ones(1,Numberof_F);% smooth tol  
%% 2系列代表筛选出的数据集
T_save2 = [];
S_save2 = [];
F_save2 = [];
dFdt_save2 = [];
dsdt_save2 = [];
dF2dt2_save2 = [];
ds2dt2_save2 = [];

%% noise related
Fsmooth_save2 = [];

Length2 = zeros(1,Numberof_F+1);
Length2lls = zeros(Numberof_F,Numberof_Tp+1);
range_save0 = zeros(Numberof_F+1,10);

%%
for i = 1:Numberof_F
    m = m_save(i,:);
    FliM = FliM_save(i,:);
    F0 = F0_save(i,:);
    F = F_save(i,:);
    s = S_save(i,:);
    Tpi = Tp(i,2:end-1);
    s_n = [0,Tpi/dt00,discreteT_length];

    %% filtering and smooth
    X = [];
    dFdt = [];
    dF2dt2 = [];

    F_smooth = [];
    F_correct = [];
    for k = 1:length(s_n)-1
        Xk = Xk_cut1 + s_n(k):s_n(k+1)-Xk_cut2; 
        Xk_cutall =  1 + s_n(k):s_n(k+1);
        X = [X,Xk]; 
        Length2lls(i,k+1) = Length2lls(i,k) + length(Xk);
        
        %% 光滑后针对完整小段
        if isempty(Xk)==0
            F_smoothk_all = csaps(t00(Xk_cutall),F(Xk_cutall),p_smooth(i),t00(Xk_cutall));
        end
        F_smoothk = F_smoothk_all(Xk_cut1:end-Xk_cut2);
        F_smooth = [F_smooth,F_smoothk];
        dFdt = [dFdt, gradient1_f(F_smoothk,dt10)];
        dF2dt2 = [dF2dt2, gradient2_f(F_smoothk,dt10)];  
        
        F_Xk = F(Xk);
        F_Xk(1:2) = F_smoothk(1:2);
        F_correct = [F_correct,F_Xk];   
        
    end

    %% update space
    T_save2 = [T_save2,t10(X)];
    S_save2 = [S_save2,s(X)];
%     F0_save2 = [F0_save2,F0(X)];
    F_save2 = [F_save2,F_correct];
    dFdt_save2 = [dFdt_save2,dFdt];
    dF2dt2_save2 = [dF2dt2_save2,dF2dt2];
    dsdt_save2 = [dsdt_save2, zeros(1,length(X))];
    ds2dt2_save2 = [ds2dt2_save2, zeros(1,length(X))];

    Fsmooth_save2 = [Fsmooth_save2,F_smooth];

    Length2(i+1) = length(T_save2);
    %% save range
    range_save0(i,:) = [min(F),max(F),min(s_origin(X)),max(s_origin(X)),...
        min(dFdt),max(dFdt),0,0,...
        min(dF2dt2),max(dF2dt2)];
    %% figures
    figure(1);   
    subplot(Fi,Fj,i);
    plot(t00,m,'-','linewidth',2); hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$m$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,FliM_save(i,:),'--');hold on
    ylabel('$FliM$','Interpreter','latex','Fontsize',18);   
    

    figure(2);   
    subplot(Fi,Fj,i);
    plot(t00,F,'-','linewidth',2); hold on
    plot(t00(X),F_smooth,'.');hold on
    plot(t00,F0,':','linewidth',2); hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$f$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,s,'--','linewidth',2);hold on
    ylabel('$s$','Interpreter','latex','Fontsize',18);

    figure(3);   
    subplot(Fi,Fj,i);
    plot(t10(X),dFdt,'.');hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$df/dt$','Interpreter','latex','Fontsize',18);

    figure(4);   
    subplot(Fi,Fj,i);
    plot(t10(X),dF2dt2,'.');hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$d^2f/dt^2$','Interpreter','latex','Fontsize',18); 
    
end
%%
range_save0(Numberof_F+1,:) = [min(range_save0(:,1)),max(range_save0(:,2)),min(range_save0(:,3)),max(range_save0(:,4)),...
    min(range_save0(:,5)),max(range_save0(:,6)),min(range_save0(:,7)),max(range_save0(:,8)),...
    min(range_save0(:,9)),max(range_save0(:,10))];

%%
save('bias2-10tpn-step1.mat','nl_snr','t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'S_para','Tp','p_smooth','m_save','FliM_save','Length2','Length2lls','T_save2',...
    'S_save','S_save2','dsdt_save2','ds2dt2_save2',...
    'F_save','F_save2','dFdt_save2','dF2dt2_save2','range_save0',...
    'Fsmooth_save2','S0_save','F0_save');