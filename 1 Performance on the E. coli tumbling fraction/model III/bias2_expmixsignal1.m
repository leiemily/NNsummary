clear;
%% smooth signals expmix
%% 条数，长度，Tp
Numberof_F = 4;
discreteT_length = 1003;
Numberof_Tp = 4;
load('bias2-10tp-mix1.mat','t_scale','dt00');

%% noise level
load('bias2-10tpn-expmix1.mat','nl_snr');
% nl_snr = unifrnd(20,22,1,Numberof_F);
snr_linear = 10.^(nl_snr/10);% 将SNR转换为线性比例

%% s setting
step_maxvalue = 2;
% rng("shuffle");
% per = unifrnd(0,0.05,1,Numberof_F);
% S_para = 0.05*unifrnd(3,5,Numberof_F,Numberof_Tp);
% Tp = [unifrnd(25,40,Numberof_F,1),unifrnd(30,45,Numberof_F,1),unifrnd(45,60,Numberof_F,1),unifrnd(55,70,Numberof_F,1)];
load('bias2-10tpn-expmix1.mat','S_para','Tp','per');
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
%% save space
S_save = zeros(Numberof_F,discreteT_length);
F_save = zeros(Numberof_F,discreteT_length);
m_save = zeros(Numberof_F,discreteT_length);
FliM_save = zeros(Numberof_F,discreteT_length);
Fs_save = zeros(Numberof_F,discreteT_length);
%% gradient save
dsdt_save = zeros(Numberof_F,discreteT_length);
dFdtorigin_save = zeros(Numberof_F,discreteT_length);
dFdt_save = zeros(Numberof_F,discreteT_length);
ds2dt2_save = zeros(Numberof_F,discreteT_length);
dF2dt2_save = zeros(Numberof_F,discreteT_length);

%% noise related
Fnoise_save = zeros(Numberof_F,discreteT_length);
Fsmooth_save = zeros(Numberof_F,discreteT_length);
dFdtsmooth_save = zeros(Numberof_F,discreteT_length);
dF2dt2smooth_save = zeros(Numberof_F,discreteT_length);


%% iteration
for i = 1:Numberof_F
   %% exp cos
    s = linearexpcom_s(S_para(i,:),Tp(i,:),t00,step_maxvalue,per(i));

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
p_smooth = 0.02*ones(1,Numberof_F);% smooth tol
p_signal = zeros(1,Numberof_F);
p_noise = zeros(1,Numberof_F);
for i = 1:Numberof_F
    F = F_save(i,:);
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
    %% smooth
    F_smooth = csaps(t00,F_noise,p_smooth(i),t00);
    
    %% update space
    Fnoise_save(i,:) = F_noise;
    Fsmooth_save(i,:) = F_smooth;
end

%% time scaling for gradient
for i = 1:Numberof_F
    s = S_save(i,:);
    F = F_save(i,:);
    F_smooth = Fsmooth_save(i,:);
    %% gradients 修改了求导方式！！！！！！！！！！！！
    dsdt = gradient1_s(s,dt10);
    dFdt = gradient1_f(F,dt10);
    dFdt_smooth = gradient1_f(F_smooth,dt10);
    
    ds2dt2 = gradient2_s(s,dt10);
    dF2dt2 = gradient2_f(F,dt10);
    dF2dt2_smooth = gradient2_f(F_smooth,dt10);

    %% update space
    dsdt_save(i,:) = dsdt;
    dFdt_save(i,:) = dFdt;
    
    ds2dt2_save(i,:) = ds2dt2;
    dF2dt2_save(i,:) = dF2dt2;

    dFdtsmooth_save(i,:) = dFdt_smooth;
    dF2dt2smooth_save(i,:) = dF2dt2_smooth;
end
%% s_smooth minus, cut
Xk_cut1 = 2;
Xk_cut2 = 1;
s_cut = Xk_cut1:length(t00)-Xk_cut2;
S_save = S_save(:,s_cut);
F_save = F_save(:,s_cut);
m_save = m_save(:,s_cut);
FliM_save = FliM_save(:,s_cut);
Fs_save = Fs_save(:,s_cut);

dsdt_save = dsdt_save(:,s_cut);
ds2dt2_save = ds2dt2_save(:,s_cut);
dFdt_save = dFdt_save(:,s_cut);
dF2dt2_save = dF2dt2_save(:,s_cut);

Fnoise_save = Fnoise_save(:,s_cut);
Fsmooth_save = Fsmooth_save(:,s_cut);
dFdtsmooth_save = dFdtsmooth_save(:,s_cut);
dF2dt2smooth_save = dF2dt2smooth_save(:,s_cut);

t00 = t00(s_cut)-t00(Xk_cut1);
t10 = t00/t_scale;
discreteT_length = length(t00);
%% save range
range_save0 = [min(F_save,[],2),max(F_save,[],2),min(S_save,[],2),max(S_save,[],2),...
    min(dFdt_save,[],2),max(dFdt_save,[],2),min(dsdt_save,[],2),max(dsdt_save,[],2),...
    min(dF2dt2_save,[],2),max(dF2dt2_save,[],2),min(ds2dt2_save,[],2),max(ds2dt2_save,[],2)];

range_save = [min(Fsmooth_save,[],2),max(Fsmooth_save,[],2),min(S_save,[],2),max(S_save,[],2),...
    min(dFdtsmooth_save,[],2),max(dFdtsmooth_save,[],2),min(dsdt_save,[],2),max(dsdt_save,[],2),...
    min(dF2dt2smooth_save,[],2),max(dF2dt2smooth_save,[],2),min(ds2dt2_save,[],2),max(ds2dt2_save,[],2)];
range_save0(Numberof_F+1,:) = [min(range_save0(:,1)),max(range_save0(:,2)),min(range_save0(:,3)),max(range_save0(:,4)),...
    min(range_save0(:,5)),max(range_save0(:,6)),min(range_save0(:,7)),max(range_save0(:,8)),...
    min(range_save0(:,9)),max(range_save0(:,10)),min(range_save0(:,11)),max(range_save0(:,12))];

range_save(Numberof_F+1,:) = [min(range_save(:,1)),max(range_save(:,2)),min(range_save(:,3)),max(range_save(:,4)),...
    min(range_save(:,5)),max(range_save(:,6)),min(range_save(:,7)),max(range_save(:,8)),...
    min(range_save(:,9)),max(range_save(:,10)),min(range_save(:,11)),max(range_save(:,12))];

%% figures
Fi = Numberof_F/2;
Fj = 2;
for i = 1:Numberof_F
%     figure(2);   
%     subplot(Fi,Fj,i);
%     plot(t00,m_save(i,:),'-','linewidth',2); hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$m$','Interpreter','latex','Fontsize',18); 
%     yyaxis right;
%     plot(t00,FliM_save(i,:),'--','linewidth',2);hold on
%     ylabel('$FliM$','Interpreter','latex','Fontsize',18); 

    figure(3);   
    subplot(Fi,Fj,i);
    plot(t00,Fnoise_save(i,:),'.'); hold on
    plot(t00,F_save(i,:),'-','linewidth',2); hold on
    plot(t00,Fsmooth_save(i,:),'-','linewidth',2); hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$f$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,S_save(i,:),'--');hold on
    ylabel('$s$','Interpreter','latex','Fontsize',18); 

    figure(4);   
    subplot(Fi,Fj,i);
    plot(t10,dFdt_save(i,:),'-','linewidth',2);hold on
    plot(t10,dFdtsmooth_save(i,:),'--','linewidth',2);hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$df/dt$','Interpreter','latex','Fontsize',18);
    yyaxis right;
    plot(t10,dsdt_save(i,:),':','linewidth',2);
    ylabel('$ds/dt$','Interpreter','latex','Fontsize',18); 

    figure(5);   
    subplot(Fi,Fj,i);
    plot(t10,dF2dt2_save(i,:),'-','linewidth',2);hold on
    plot(t10,dF2dt2smooth_save(i,:),'--','linewidth',2);hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$d^2f/dt^2$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t10,ds2dt2_save(i,:),':','linewidth',2);
    ylabel('$d^2s/dt^2$','Interpreter','latex','Fontsize',18); 
end

%% % 验证精度
error_noisef = zeros(1,Numberof_F);
error_noisedf = zeros(1,Numberof_F);
error_noisedf2 = zeros(1,Numberof_F);
for i = 1:Numberof_F
    error_noisedf2(i) = norm(dF2dt2_save(i,:) - dF2dt2smooth_save(i,:),2)/norm(dF2dt2_save(i,:),2);
    error_noisedf(i) = norm(dFdt_save(i,:) - dFdtsmooth_save(i,:),2)/norm(dFdt_save(i,:),2);
    error_noisef(i) = mean((F_save(i,:) - Fsmooth_save(i,:)).^2./F_save(i,:).^2);
end

range_err = [error_noisedf2;error_noisedf;error_noisef]';

%%
save('bias2-10tp-expmix1.mat','nl_snr','per','t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'S_para','Tp','p_smooth','m_save','FliM_save','Fs_save',...
    'S_save','dsdt_save','ds2dt2_save',...
    'F_save','dFdt_save','dF2dt2_save',...
    'Fnoise_save','Fsmooth_save','dFdtsmooth_save','dF2dt2smooth_save',...
    'range_save0','range_save','range_err');