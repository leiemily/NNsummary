clear;
%% step signals
%% 条数，长度，Tp
Numberof_F = 16;
discreteT_length = 1071;
Numberof_Tp = 10;
% dt00 = 0.2;
load('frac2-1tp-mix1.mat','t_scale','dt00');

%% s setting
% load('frac2-1tp-step1.mat','S_para','Tp');

step_maxvalue = 2;
step_minlen = (discreteT_length-1)*dt00/15;

rng("shuffle");
S_para = zeros(Numberof_F,Numberof_Tp);
for i = 1:Numberof_F
    Si_diff = step_maxvalue*ones(1,Numberof_Tp);
    while max(abs(Si_diff)) > step_maxvalue/2.2
        Si = unifrnd(0,step_maxvalue,1, Numberof_Tp);
        Si_diff = diff(Si);
    end
    S_para(i,:) = Si;
end

% S_para = unifrnd(0,step_maxvalue,Numberof_F, Numberof_Tp-1);
% S_para = [zeros(Numberof_F,1),S_para];

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
% load('frac2-10tp-step1.mat','S_para','Tp');

%% figures
Fi = Numberof_F/2;
Fj = 2;
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
t10 = t00/t_scale;
dt10 = dt00/t_scale;

%% cut
Xk_cut1 = 7;
Xk_cut2 = 1;
%% save space
S_save = zeros(Numberof_F,discreteT_length);
F_save = zeros(Numberof_F,discreteT_length);
m_save = zeros(Numberof_F,discreteT_length);
Yp_save = zeros(Numberof_F,discreteT_length);
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
    F = cellmovemodel_3(Yp);
    %% update space
    S_save(i,:) = s;
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
    
    s_n = [0,find(diff(s)),discreteT_length];
   
    
    %% filtering and smooth
    X = [];
    dFdt = [];
    dF2dt2 = [];

    for k = 1:length(s_n)-1
        Xk = Xk_cut1 + s_n(k):s_n(k+1)-Xk_cut2; 
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
    dsdt_save2 = [dsdt_save2, zeros(1,length(X))];
    ds2dt2_save2 = [ds2dt2_save2, zeros(1,length(X))];

    Length2(i+1) = length(T_save2);
    %% save range
    range_save0(i,:) = [min(F),max(F),min(s(X)),max(s(X)),...
        min(dFdt),max(dFdt),0,0,...
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
save('frac2-1tp-step1.mat','t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'S_para','Tp','m_save','Yp_save','Length2','Length2lls','T_save2',...
    'S_save','S_save2','dsdt_save2','ds2dt2_save2',...
    'F_save','F_save2','dFdt_save2','dF2dt2_save2',...
    'range_save0');