clear;
%% step signals
load('tumbling_fraction_step1all');

%% 其中8条
Fnoise_save = Fnoise_save([1,6,7,8],:);
S_save = S_save([1,6,7,8],:);

%% 截断点 
part_length = 401;
s_cut = 1:part_length;%length(t);
part_start = [62,470,993,1400,1950,2615,3020;...
    10,811,1187,1891,2295,2690,3105;...
    1,420,1000,2041,2420,2830,3125;...
    396,821,1206,1611,2071,2594,3064];
% part_start = [62,479,993,1400,1950,2615,3020;...
%     10,811,1187,1891,2295,2700,3105;...
%     1,420,1000,2011,2420,2830,3125;...
%     396,801,1206,1611,2071,2594,3044];
f_save = zeros(size(part_start,1)*size(part_start,2),part_length);
s_save = zeros(size(part_start,1)*size(part_start,2),part_length);
for i = 1:size(part_start,1)
    for j = 1:size(part_start,2)
        f_save((i-1)*size(part_start,2)+j,:) = Fnoise_save(i,part_start(i,j):part_start(i,j)+part_length-1);
        s_save((i-1)*size(part_start,2)+j,:) = S_save(i,part_start(i,j):part_start(i,j)+part_length-1);
    end
end
S_save = s_save;
Fnoise_save = f_save;

% Numberof_F = size(Fnoise_save,1);
% Numberof_Tp = 20;

%% rescaling
t_scale = 200;
dt00 = 1;
t00 = t(s_cut) - 0.9;%改成1:...
t10 = t00/t_scale;
dt10 = dt00/t_scale;
discreteT_length = length(t00);

%% 记录条数
Numberof_F = size(Fnoise_save,1);
Numberof_Tp = 5;

%% figures
Fi = 4;
Fj = Numberof_F/Fi;

%% smooth tol
p_smooth = 5e-5*ones(1,Numberof_F);
% length_cut = 100;
Xk_cut1 = 4;
Xk_cut2 = 1;

%% save space
T_save2 = [];
S_save2 = [];
dsdt_save2 = [];
ds2dt2_save2 = [];

%% noise save
Fnoise_save2 = [];
Fsmooth_save2 = [];
dFdtsmooth_save2 = [];
dF2dt2smooth_save2 = [];

Length2 = zeros(1,Numberof_F+1);
Length2lls = zeros(Numberof_F,Numberof_Tp);
range_save = zeros(Numberof_F+1,10);
%% iter
for i = 1:Numberof_F
    s = S_save(i,:);
    F_noise = Fnoise_save(i,:);

    s_n = [0,find(diff(s)),discreteT_length];
    %% filtering and smooth
    X = [];
    F_smooth = [];
    dFdt_smooth = [];
    dF2dt2_smooth = [];
    nn = 1;
    for k = 1:length(s_n)-1
        Xk = Xk_cut1 + s_n(k):s_n(k+1) - Xk_cut2;
        X = [X,Xk]; 
        Length2lls(i,nn+1) = Length2lls(i,nn) + length(Xk);
        nn = nn+1;

        Xk_cutall =  1 + s_n(k):s_n(k+1);
        F_smoothk = csaps(t00(Xk_cutall),F_noise(Xk_cutall),p_smooth(i),t00(Xk_cutall));
        dFdt_smoothk = gradient1_f(F_smoothk(Xk_cut1:end-Xk_cut2),dt10);
        dF2dt2_smoothk = gradient2_f(F_smoothk(Xk_cut1:end-Xk_cut2),dt10);

        F_smooth = [F_smooth,F_smoothk(Xk_cut1:end-Xk_cut2)];
        dFdt_smooth = [dFdt_smooth, dFdt_smoothk];
        dF2dt2_smooth = [dF2dt2_smooth, dF2dt2_smoothk];   
            
    end
    %% update space
    T_save2 = [T_save2,t10(X)];
    S_save2 = [S_save2,s(X)];
    dsdt_save2 = [dsdt_save2, zeros(1,length(X))];
    ds2dt2_save2 = [ds2dt2_save2, zeros(1,length(X))];

    Fnoise_save2 = [Fnoise_save2,F_noise(X)];
    Fsmooth_save2 = [Fsmooth_save2,F_smooth];
    dFdtsmooth_save2 = [dFdtsmooth_save2, dFdt_smooth];
    dF2dt2smooth_save2 = [dF2dt2smooth_save2, dF2dt2_smooth];

    Length2(i+1) = length(T_save2);
    %% save range
    range_save(i,:) = [min(F_smooth),max(F_smooth),min(s(X)),max(s(X)),...
        min(dFdt_smooth),max(dFdt_smooth),0,0,...
        min(dF2dt2_smooth),max(dF2dt2_smooth)];
    %% figures
    figure(2);   
    subplot(Fi,Fj,i);
    plot(t00,F_noise,'.'); hold on
    plot(t00(X),F_smooth,'.');hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$f$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,s,'--','linewidth',2);hold on
    ylabel('$s$','Interpreter','latex','Fontsize',18);

    figure(3);   
    subplot(Fi,Fj,i);
    plot(t10(X),dFdt_smooth,'.');hold on
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$df/dt$','Interpreter','latex','Fontsize',18);

    figure(4);   
    subplot(Fi,Fj,i);
    plot(t10(X),dF2dt2_smooth,'.');
    xlabel('$t$','Interpreter','latex','Fontsize',18);
    ylabel('$d^2f/dt^2$','Interpreter','latex','Fontsize',18); 
end
range_save(Numberof_F+1,:) = [min(range_save(:,1)),max(range_save(:,2)),min(range_save(:,3)),max(range_save(:,4)),...
    min(range_save(:,5)),max(range_save(:,6)),min(range_save(:,7)),max(range_save(:,8)),...
    min(range_save(:,9)),max(range_save(:,10))];

%%
save('exp1-0tp-step1.mat','t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'p_smooth','Length2','Length2lls','T_save2',...
    'S_save','S_save2','dsdt_save2','ds2dt2_save2',...
    'Fnoise_save','Fsmooth_save2','dFdtsmooth_save2','dF2dt2smooth_save2',...
    'range_save');