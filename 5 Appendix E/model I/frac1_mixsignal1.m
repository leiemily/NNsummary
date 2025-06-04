clear;
%% smooth signals mix
%% 条数，长度，Tp
Numberof_F = 10;
discreteT_length = 1003;
Numberof_Tp = 2;
dt00 = 0.2;

%% noise level
nl_snr = unifrnd(20,22,1,Numberof_F);
snr_linear = 10.^(nl_snr/10);% 将SNR转换为线性比例


%% s setting
% rng("shuffle");
% S_para = 0.05*unifrnd(2,5,Numberof_F,Numberof_Tp);
% Tp = 10*[unifrnd(1,2,Numberof_F,1),unifrnd(2,3,Numberof_F,1)];

load('frac1-10tpn-mix1.mat','S_para','Tp');

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

%% save space
S_save = zeros(Numberof_F,discreteT_length);%储存加噪s
F_save = zeros(Numberof_F,discreteT_length);%储存加噪s下的f
S0_save = zeros(Numberof_F,discreteT_length);%储存无噪s
F0_save = zeros(Numberof_F,discreteT_length);%储存无噪s下的f

m_save = zeros(Numberof_F,discreteT_length);
Yp_save = zeros(Numberof_F,discreteT_length);
%% gradient save
dsdt_save = zeros(Numberof_F,discreteT_length);
dFdt_save = zeros(Numberof_F,discreteT_length);
ds2dt2_save = zeros(Numberof_F,discreteT_length);
dF2dt2_save = zeros(Numberof_F,discreteT_length);

%% noise related 确认光滑化效果
Fsmooth_save = zeros(Numberof_F,discreteT_length);
Ssmooth_save = zeros(Numberof_F,discreteT_length);

%% iteration
%% noiseless
for i = 1:Numberof_F
    %% sin linear combination signal s(t)(changed)
    s = linearcom_s(S_para(i,:),Tp(i,:),t00,Numberof_Tp);

    %% save space
    m = zeros(1,discreteT_length);
    A = zeros(1,discreteT_length);
    E = zeros(1,discreteT_length);

    R = zeros(1,discreteT_length);
    Bp = zeros(1,discreteT_length);
    Yp = zeros(1,discreteT_length);

    %% initial space
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
    S0_save(i,:) = s;
    F0_save(i,:) = F;
    
end

%% noisy
p_signal = zeros(1,Numberof_F);
p_noise = zeros(1,Numberof_F);

for i = 1:Numberof_F
    %% sin linear combination signal s(t)(changed)
    s_origin = S0_save(i,:);
    %% add  noise
    p_signal(i) = mean(s_origin.^2);
    p_noise(i) = p_signal(i) / snr_linear(i);% 计算噪声功率
    
    sigma = sqrt(p_noise(i));
    
    rng('default');
    seed_noise = rng;
    rng(seed_noise);
    noise_s = normrnd(0,sigma,size(s_origin));
%     noise_s(1:3) = 0;%初始不含噪
    s = s_origin + noise_s;
    s_overnn = find(s<0);
    s(s_overnn) = s_origin(s_overnn);

    %% save space
    m = zeros(1,discreteT_length);
    A = zeros(1,discreteT_length);
    E = zeros(1,discreteT_length);

    R = zeros(1,discreteT_length);
    Bp = zeros(1,discreteT_length);
    Yp = zeros(1,discreteT_length);

    %% initial space
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
    S_save(i,:) = s;
    F_save(i,:) = F;
    
    m_save(i,:) = m;
    Yp_save(i,:) = Yp;
end

%% noise and smooth noise and smooth  noise and smooth  noise and smooth  noise and smooth  
p_smooth = [0.3,0.3,0.3,0.3,0.3,...
            0.3,0.3,0.4,0.3,0.4];

%% time scaling for gradient
% t_scale = round(1/max(dFdtorigin_save,[],'all'),1);
t_scale = 6.4;
dt10 = dt00/t_scale;

for i = 1:Numberof_F
    F = F_save(i,:);
    s = S_save(i,:);
    %% smooth
    F_smooth = csaps(t00,F,p_smooth(i),t00);
    S_smooth = csaps(t00,s,p_smooth(i),t00);
    Ssmooth_save(i,:) = S_smooth;
    Fsmooth_save(i,:) = F_smooth;
    

    %% gradients 修改了求导方式！！！！！！！！！！！！
    dsdt = gradient1_s(S_smooth,dt10);
    dFdt = gradient1_f(F_smooth,dt10);
    
    ds2dt2 = gradient2_s(S_smooth,dt10);
    dF2dt2 = gradient2_f(F_smooth,dt10);
    
    dsdt_save(i,:) = dsdt;
    dFdt_save(i,:) = dFdt;

    ds2dt2_save(i,:) = ds2dt2;
    dF2dt2_save(i,:) = dF2dt2;
end


%% s_smooth minus, cut
Xk_cut1 = 2;
Xk_cut2 = 1;
s_cut = Xk_cut1:length(t00)-Xk_cut2;
m_save = m_save(:,s_cut);
Yp_save = Yp_save(:,s_cut);

Ssmooth_save = Ssmooth_save(:,s_cut);
Fsmooth_save = Fsmooth_save(:,s_cut);
S0_save = S0_save(:,s_cut);
F0_save = F0_save(:,s_cut);
S_save = S_save(:,s_cut);
F_save = F_save(:,s_cut);
F_save(:,1:2) = Fsmooth_save(:,1:2);%f初始精度高才行

dsdt_save = dsdt_save(:,s_cut);
ds2dt2_save = ds2dt2_save(:,s_cut);
dFdt_save = dFdt_save(:,s_cut);
dF2dt2_save = dF2dt2_save(:,s_cut);

t00 = t00(s_cut)-t00(Xk_cut1);
t10 = t00/t_scale;
discreteT_length = length(t00);
%% save range
range_save0 = [min(F_save,[],2),max(F_save,[],2),min(S_save,[],2),max(S_save,[],2),...
    min(dFdt_save,[],2),max(dFdt_save,[],2),min(dsdt_save,[],2),max(dsdt_save,[],2),...
    min(dF2dt2_save,[],2),max(dF2dt2_save,[],2),min(ds2dt2_save,[],2),max(ds2dt2_save,[],2)];
range_save0(Numberof_F+1,:) = [min(range_save0(:,1)),max(range_save0(:,2)),min(range_save0(:,3)),max(range_save0(:,4)),...
    min(range_save0(:,5)),max(range_save0(:,6)),min(range_save0(:,7)),max(range_save0(:,8)),...
    min(range_save0(:,9)),max(range_save0(:,10)),min(range_save0(:,11)),max(range_save0(:,12))];


%% % figures
Fi = Numberof_F/2;
Fj = 2;
for i = 1:Numberof_F
    figure(2);   
    subplot(Fi,Fj,i);
    plot(t00,m_save(i,:),'-','linewidth',2); hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$m$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,Yp_save(i,:),'--');hold on
    ylabel('$Yp$','Interpreter','latex','Fontsize',18); 
    
    
    figure(3);   
    subplot(Fi,Fj,i);
    plot(t00,F_save(i,:),'.'); hold on
    plot(t00,Fsmooth_save(i,:),'-'); hold on
    plot(t00,F0_save(i,:),':'); hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$f$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,S_save(i,:),'--');hold on
    plot(t00,Ssmooth_save(i,:),'-'); hold on
    plot(t00,S0_save(i,:),':');hold on
    ylabel('$s$','Interpreter','latex','Fontsize',18); 

    figure(4);   
    subplot(Fi,Fj,i);
    plot(t10,dFdt_save(i,:),'-','linewidth',2);hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$df/dt$','Interpreter','latex','Fontsize',18);
    yyaxis right;
    plot(t10,dsdt_save(i,:),':','linewidth',2);
    ylabel('$ds/dt$','Interpreter','latex','Fontsize',18); 

    figure(5);   
    subplot(Fi,Fj,i);
    plot(t10,dF2dt2_save(i,:),'-','linewidth',2);hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$d^2f/dt^2$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t10,ds2dt2_save(i,:),':','linewidth',2);
    ylabel('$d^2s/dt^2$','Interpreter','latex','Fontsize',18); 
end

%%
save('frac1-10tpn-mix1.mat','nl_snr','t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'S_para','Tp','p_smooth','m_save','Yp_save',...
    'S_save','dsdt_save','ds2dt2_save',...
    'F_save','dFdt_save','dF2dt2_save','range_save0',...
    'Ssmooth_save','Fsmooth_save','S0_save','F0_save');
%     'Fnoise_save','Fsmooth_save','dFdtsmooth_save','dF2dt2smooth_save',...
%     'range_save0','range_save','range_err');