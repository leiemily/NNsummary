clear;
%% smooth signals mix
%% 条数，长度，Tp
Numberof_F = 10;
discreteT_length = 1003;
Numberof_Tp = 2;
% dt00 = 0.2;
load('frac2-10tp-mix1.mat','dt00');

%% noise level
load('frac2-10tp-mix1.mat','nl_snr');
% nl_snr = unifrnd(20,22,1,Numberof_F);
snr_linear = 10.^(nl_snr/10);% 将SNR转换为线性比例


%% s setting
load('frac2-10tp-mix1','S_para','Tp','S_save1uncut');

% %% 9tp
% stretch_para = 9;
% p_smooth = 0.2*ones(1,Numberof_F);

% %% 8tp
% stretch_para = 8;
% p_smooth = 0.3*ones(1,Numberof_F);

% %% 7tp
% stretch_para = 7;
% p_smooth = 0.4*ones(1,Numberof_F);

% %% 6tp
% stretch_para = 6;
% p_smooth = 0.45*ones(1,Numberof_F);

% %% 5tp
% stretch_para = 5;
% p_smooth = 0.5*ones(1,Numberof_F);

%% 4tp
stretch_para = 4;
p_smooth = 0.75*ones(1,Numberof_F);

%%
dt00 = stretch_para*dt00/10;
Tp = stretch_para*Tp/10;
S_para = stretch_para*S_para/10;

%% para
kR = 3.82e-2;
a0 = 0.5;
m0 = 1;
Nr = 6;
alpha0 = 1.7;
KI = 18.2;
KA = 3000;

ka = 10;%未找到文献数据
tauz = 0.5;

%% rescaling
t00 = 0:dt00:(discreteT_length-1)*dt00;

%% save space
S_save = S_save1uncut;
F_save = zeros(Numberof_F,discreteT_length);
m_save = zeros(Numberof_F,discreteT_length);
Yp_save = zeros(Numberof_F,discreteT_length);

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
    %% sin linear combination signal s(t)(changed)
    s = S_save(i,:);

    %% save space
    m = zeros(1,discreteT_length);
    Yp = zeros(1,discreteT_length);
    A = zeros(1,discreteT_length);
    E = zeros(1,discreteT_length);

    %% initial space
    x0 = 2;
    c = s(1);
    fun = @(x)paramfun_m(x,c);
    m_initial = fsolve(fun,x0);

    m(1) = m_initial;
    E(1) = - alpha0 *(m(1) - m0) + log((1 + s(1)/KI)/(1 + s(1)/KA));
    A(1) = (1 + exp(Nr * E(1)))^(-1);
    Yp(1) = ka*tauz*A(1);

    %% iteration
    for j = 1:discreteT_length-1
        m(j+1) = m(j) + dt00 * kR * (1 - A(j)/a0);
        Yp(j+1) = Yp(j) + dt00 * (ka * A(j) - Yp(j)/tauz);
        E(j+1) = - alpha0 *(m(j+1) - m0) + log((1 + s(j+1)/KI)/(1 + s(j+1)/KA));
        A(j+1) = (1 + exp(Nr * E(j+1)))^(-1);
    end

    %% calculate F
    F = cellmovemodel_3(Yp);
    %% update space
    S_save(i,:) = s;
    F_save(i,:) = F;
    m_save(i,:) = m;
    Yp_save(i,:) = Yp;

    dFdt = gradient1_f(F,dt00);
    dFdtorigin_save(i,:) = dFdt;
end

%% noise and smooth
% p_smooth = 0.1*ones(1,Numberof_F);
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
t_scale = round(1/max(dFdtorigin_save,[],'all'),1);
dt10 = dt00/t_scale;
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
% S_save1uncut = S_save;
S_save = S_save(:,s_cut);
F_save = F_save(:,s_cut);
m_save = m_save(:,s_cut);
Yp_save = Yp_save(:,s_cut);

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

%% % 验证精度
error_noisef = zeros(1,Numberof_F);
error_noisedf = zeros(1,Numberof_F);
error_noisedf2 = zeros(1,Numberof_F);
for i = 1:Numberof_F
    error_noisedf2(i) = norm(dF2dt2_save(i,:) - dF2dt2smooth_save(i,:),2)/norm(dF2dt2_save(i,:),2);
    error_noisedf(i) = norm(dFdt_save(i,:) - dFdtsmooth_save(i,:),2)/norm(dFdt_save(i,:),2);
    error_noisef(i) = sum((F_save(i,:) - Fsmooth_save(i,:)).^2./F_save(i,:).^2)/discreteT_length;
end

range_err = [error_noisedf2;error_noisedf;error_noisef]';

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
    plot(t00,Yp_save(i,:),'--','linewidth',2);hold on
    ylabel('$Yp$','Interpreter','latex','Fontsize',18); 

    figure(3);   
    subplot(Fi,Fj,i);
    plot(t00,Fnoise_save(i,:),'.'); hold on
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

%%
save(['frac2-',num2str(stretch_para),'tp-mix1'],'t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'S_para','Tp','p_smooth','m_save','Yp_save','S_save1uncut',...
    'S_save','dsdt_save','ds2dt2_save',...
    'F_save','dFdt_save','dF2dt2_save',...
    'Fnoise_save','Fsmooth_save','dFdtsmooth_save','dF2dt2smooth_save',...
    'range_save0','range_save','range_err');