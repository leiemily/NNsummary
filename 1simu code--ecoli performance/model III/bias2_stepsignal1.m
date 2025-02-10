clear;
%% step signals
%% 条数，长度，Tp
Numberof_F = 10;
discreteT_length = 1002;
Numberof_Tp = 10;
% dt00 = 0.2;
load('bias2-10tp-mix1.mat','t_scale','dt00');
%% noise level
% nl_snr = unifrnd(20,22,1,Numberof_F);
nl_snr = 22*ones(1,Numberof_F);
snr_linear = 10.^(nl_snr/10);% 将SNR转换为线性比例

%% s setting
step_maxvalue = 2;
step_minlen = (discreteT_length-1)*dt00/15;

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
dt10 = dt00/t_scale;
t10 = t00/t_scale;
%% cut
Xk_cut1 = 2;
Xk_cut2 = 1;

%% save space
S_save = zeros(Numberof_F,discreteT_length);
F_save = zeros(Numberof_F,discreteT_length);
m_save = zeros(Numberof_F,discreteT_length);
FliM_save = zeros(Numberof_F,discreteT_length);
Fs_save = zeros(Numberof_F,discreteT_length);
Fnoise_save = zeros(Numberof_F,discreteT_length);

%% iteration
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
    Fs = zeros(1,discreteT_length);
    
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
    Fs(1) = (((Nr*bb*exp(Nr*E(1))*(FliM(1) - aa)*(1/(KI*(s(1)/KA + 1)) - (s(1)/KI + 1)/(KA*(s(1)/KA + 1)^2))*(1/(K_A*(exp(Nr*E(1)) + 1)) + 1)^(bb*(FliM(1) - aa) - 1)*(s(1)/KA + 1))/(K_A*(exp(Nr*E(1)) + 1)^2*(s(1)/KI + 1)) + (Nr*P*bb*exp(Nr*E(1))*(FliM(1) - aa)*(1/(C*K_A*(exp(Nr*(log((s(1)/KI + 1)/(s(1)/KA + 1))- alpha0*(m(1) - m0))) + 1)) + 1)^(bb*(FliM(1) - aa) - 1)*(1/(KI*(s(1)/KA + 1)) - (s(1)/KI + 1)/(KA*(s(1)/KA + 1)^2))*(s(1)/KA + 1))/(C*K_A*(exp(Nr*E(1)) + 1)^2*(s(1)/KI + 1)))*(1/(K_A*(exp(Nr*E(1)) + 1)) + 1)^(bb*(FliM(1) - aa)))/((1/(K_A*(exp(Nr*E(1)) + 1)) + 1)^(bb*(FliM(1) - aa))...
    + P*(1/(C*K_A*(exp(Nr*E(1)) + 1)) + 1)^(bb*(FliM(1) - aa)))^2 - (Nr*bb*exp(Nr*E(1))*(FliM(1) - aa)*(1/(KI*(s(1)/KA + 1)) - (s(1)/KI + 1)/(KA*(s(1)/KA + 1)^2))*(1/(K_A*(exp(Nr*E(1)) + 1)) + 1)^(bb*(FliM(1) - aa) - 1)*(s(1)/KA + 1))/(K_A*(exp(Nr*E(1)) + 1)^2*((1/(K_A*(exp(Nr*E(1)) + 1)) + 1)^(bb*(FliM(1) - aa))...
    + P*(1/(C*K_A*(exp(Nr*E(1)) + 1)) + 1)^(bb*(FliM(1) - aa)))*(s(1)/KI + 1));

    %% iteration
    for j = 1:discreteT_length-1
        m(j+1)  = m(j) + dt00 * (V_R*(1 - A(j))- V_B*A(j));
        FliM(j+1) = FliM(j) + dt00 * (V_on * (FliM_max - FliM(j)) - V_off * A(j)* FliM(j));

        E(j+1) = - alpha0 *(m(j+1) - m0) + log((1 + s(j+1)/KI)/(1 + s(j+1)/KA));
        A(j+1) = (1 + exp (Nr * E(j+1))).^(-1);
        Fs(j+1) = (((Nr*bb*exp(Nr*E(j+1))*(FliM(j+1) - aa)*(1/(KI*(s(j+1)/KA + 1)) - (s(j+1)/KI + 1)/(KA*(s(j+1)/KA + 1)^2))*(1/(K_A*(exp(Nr*E(j+1)) + 1)) + 1)^(bb*(FliM(j+1) - aa) - 1)*(s(j+1)/KA + 1))/(K_A*(exp(Nr*E(j+1)) + 1)^2*(s(j+1)/KI + 1)) + (Nr*P*bb*exp(Nr*E(j+1))*(FliM(j+1) - aa)*(1/(C*K_A*(exp(Nr*(log((s(j+1)/KI + 1)/(s(j+1)/KA + 1))- alpha0*(m(j+1) - m0))) + 1)) + 1)^(bb*(FliM(j+1) - aa) - 1)*(1/(KI*(s(j+1)/KA + 1)) - (s(j+1)/KI + 1)/(KA*(s(j+1)/KA + 1)^2))*(s(j+1)/KA + 1))/(C*K_A*(exp(Nr*E(j+1)) + 1)^2*(s(j+1)/KI + 1)))*(1/(K_A*(exp(Nr*E(j+1)) + 1)) + 1)^(bb*(FliM(j+1) - aa)))/((1/(K_A*(exp(Nr*E(j+1)) + 1)) + 1)^(bb*(FliM(j+1) - aa))...
        + P*(1/(C*K_A*(exp(Nr*E(j+1)) + 1)) + 1)^(bb*(FliM(j+1) - aa)))^2 - (Nr*bb*exp(Nr*E(j+1))*(FliM(j+1) - aa)*(1/(KI*(s(j+1)/KA + 1)) - (s(j+1)/KI + 1)/(KA*(s(j+1)/KA + 1)^2))*(1/(K_A*(exp(Nr*E(j+1)) + 1)) + 1)^(bb*(FliM(j+1) - aa) - 1)*(s(j+1)/KA + 1))/(K_A*(exp(Nr*E(j+1)) + 1)^2*((1/(K_A*(exp(Nr*E(j+1)) + 1)) + 1)^(bb*(FliM(j+1) - aa))...
        + P*(1/(C*K_A*(exp(Nr*E(j+1)) + 1)) + 1)^(bb*(FliM(j+1) - aa)))*(s(j+1)/KI + 1));
    end

    %% calculate F
    F = (1 + A/K_A).^((FliM-aa)*bb) ./ ((1 + A/K_A).^((FliM-aa)*bb) + P * (1 + A/(K_A*C)).^((FliM-aa)*bb));
    %% update space
    S_save(i,:) = s;
    F_save(i,:) = F;
    m_save(i,:) = m;
    FliM_save(i,:) = FliM;
    Fs_save(i,:) = Fs;
end

%% noise and smooth noise and smooth noise and smooth noise and smooth noise and smooth noise and smooth noise and smooth noise and smooth
p_smooth = 0.005*ones(1,Numberof_F);% smooth tol
p_signal = zeros(1,Numberof_F);
p_noise = zeros(1,Numberof_F);
%% 2系列代表筛选出的数据集
T_save2 = [];
S_save2 = [];
F_save2 = [];
dFdt_save2 = [];
dF2dt2_save2 = [];
dsdt_save2 = [];
ds2dt2_save2 = [];

%% noise related
Fnoise_save2 = [];
Fsmooth_save2 = [];
dFdtsmooth_save2 = [];
dF2dt2smooth_save2 = [];

Length2 = zeros(1,Numberof_F+1);
Length2lls = zeros(Numberof_F,2*Numberof_Tp);
range_save0 = zeros(Numberof_F+1,10);
range_save = zeros(Numberof_F+1,10);

error_noisef = zeros(1,Numberof_F);
error_noisedf = zeros(1,Numberof_F);
error_noisedf2 = zeros(1,Numberof_F);
%%
for i = 1:Numberof_F
    m = m_save(i,:);
    FliM = FliM_save(i,:);
    F = F_save(i,:);
    s = S_save(i,:);
    s_n = [0,find(diff(s)),discreteT_length];
    %% add  noise
    p_signal(i) = mean(F.^2);
    p_noise(i) = p_signal(i) / snr_linear(i);% 计算噪声功率
    
    sigma = sqrt(p_noise(i));
    
    rng('default');
    seed_noise = rng;
    rng(seed_noise);
    noise_F = normrnd(0,sigma,size(F));
    F_noise = F + noise_F;
    %% constraint F_noise
    F_overnn = find(F_noise>1 | F_noise<0);
    F_noise(F_overnn) = F(F_overnn);
   
    %% filtering and smooth
    X = [];
    dFdt = [];
    dF2dt2 = [];

    F_smooth = [];
    dFdt_smooth = [];
    dF2dt2_smooth = [];
    for k = 1:length(s_n)-1
        Xk = Xk_cut1 + s_n(k):s_n(k+1)-Xk_cut2; 
        Xk_cutall =  1 + s_n(k):s_n(k+1);
        X = [X,Xk]; 
        Length2lls(i,k+1) = Length2lls(i,k) + length(Xk);
        
        %% 光滑后针对完整小段
        if isempty(Xk)==0
            F_smoothk_all = csaps(t00(Xk_cutall),F_noise(Xk_cutall),p_smooth(i),t00(Xk_cutall));
        end
        F_smoothk = F_smoothk_all(Xk_cut1:end-Xk_cut2);
        F_smooth = [F_smooth,F_smoothk];
        dFdt = [dFdt, gradient1_f(F(Xk),dt10)];
        dFdt_smooth = [dFdt_smooth, gradient1_f(F_smoothk,dt10)];
        dF2dt2 = [dF2dt2, gradient2_f(F(Xk),dt10)];
        dF2dt2_smooth = [dF2dt2_smooth, gradient2_f(F_smoothk,dt10)];     
        
    end

    error_noisef(i) = sum((F(X) - F_smooth).^2/F(X).^2)/length(X);
    error_noisedf(i) = norm(dFdt - dFdt_smooth,2)/norm(dFdt,2);
    error_noisedf2(i) = norm(dF2dt2 - dF2dt2_smooth,2)/norm(dF2dt2,2);
    
    %% update space
    Fnoise_save(i,:) = F_noise;
    
    T_save2 = [T_save2,t10(X)];
    S_save2 = [S_save2,s(X)];
    F_save2 = [F_save2,F(X)];
    dFdt_save2 = [dFdt_save2,dFdt];
    dF2dt2_save2 = [dF2dt2_save2,dF2dt2];
    dsdt_save2 = [dsdt_save2, zeros(1,length(X))];
    ds2dt2_save2 = [ds2dt2_save2, zeros(1,length(X))];

    Fnoise_save2 = [Fnoise_save2,F_noise(X)];
    Fsmooth_save2 = [Fsmooth_save2,F_smooth];
    dFdtsmooth_save2 = [dFdtsmooth_save2, dFdt_smooth];
    dF2dt2smooth_save2 = [dF2dt2smooth_save2, dF2dt2_smooth];

    Length2(i+1) = length(T_save2);
   %% save range
    range_save0(i,:) = [min(F),max(F),min(s(X)),max(s(X)),...
        min(dFdt),max(dFdt),0,0,...
        min(dF2dt2),max(dF2dt2)];
    
    range_save(i,:) = [min(F_smooth),max(F_smooth),min(s(X)),max(s(X)),...
        min(dFdt_smooth),max(dFdt_smooth),0,0,...
        min(dF2dt2_smooth),max(dF2dt2_smooth)];
   %% figures
%     figure(1);   
%     subplot(Fi,Fj,i);
%     plot(t00,m,'-','linewidth',2); hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$m$','Interpreter','latex','Fontsize',18); 
%     yyaxis right;
%     plot(t00,FliM,'--','linewidth',2);hold on
%     ylabel('$FliM$','Interpreter','latex','Fontsize',18); 
    
    figure(2);   
    subplot(Fi,Fj,i);
    plot(t00,F_noise,'.'); hold on
    plot(t00,F,'-','linewidth',2); hold on
    plot(t00(X),F_smooth,'.');hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$f$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,s,'--','linewidth',2);hold on
    ylabel('$s$','Interpreter','latex','Fontsize',18);

    figure(3);   
    subplot(Fi,Fj,i);
    plot(t10(X),dFdt,'.');hold on
    plot(t10(X),dFdt_smooth,'.');hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$df/dt$','Interpreter','latex','Fontsize',18);

    figure(5);   
    subplot(Fi,Fj,i);
    plot(t10(X),dF2dt2,'.');hold on
    plot(t10(X),dF2dt2_smooth,'.');
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$d^2f/dt^2$','Interpreter','latex','Fontsize',18); 
end
range_save0(Numberof_F+1,:) = [min(range_save0(:,1)),max(range_save0(:,2)),min(range_save0(:,3)),max(range_save0(:,4)),...
    min(range_save0(:,5)),max(range_save0(:,6)),min(range_save0(:,7)),max(range_save0(:,8)),...
    min(range_save0(:,9)),max(range_save0(:,10))];

range_save(Numberof_F+1,:) = [min(range_save(:,1)),max(range_save(:,2)),min(range_save(:,3)),max(range_save(:,4)),...
    min(range_save(:,5)),max(range_save(:,6)),min(range_save(:,7)),max(range_save(:,8)),...
    min(range_save(:,9)),max(range_save(:,10))];

range_err = [error_noisedf2;error_noisedf;error_noisef]';

%%
save('bias2-10tp-step1.mat','nl_snr','t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'p_smooth','m_save','FliM_save','Length2','Length2lls','T_save2',...
    'S_para','Tp','S_save','S_save2','dsdt_save2','ds2dt2_save2',...
    'F_save','F_save2','dFdt_save2','dF2dt2_save2',...
    'Fnoise_save','Fsmooth_save2','dFdtsmooth_save2','dF2dt2smooth_save2',...
    'range_save0','range_save','range_err');