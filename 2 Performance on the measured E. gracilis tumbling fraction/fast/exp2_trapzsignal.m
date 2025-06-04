clear;
load('tumbling_fraction_trapzall');

%% trapz signals
%% 截断点 
part_length = 404;
s_cut = 1:part_length;%length(t);

part_start = [184,630,1429,2256;380,997,1455,2400;325,776,1422,2800];
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



%% time 
t_scale = 100;
dt00 = 1;
t00 = t(s_cut)-0.9;
dt10 = dt00/t_scale;
t10 = t00/t_scale;
discreteT_length = length(t10);

%% 记录条数
Numberof_F = size(Fnoise_save,1);
Numberof_Tp = 8;

%% 调整顺序
A_trapz = 1:Numberof_F;
C_trapz = [3,5,10];
B_trapz = setdiff(A_trapz,C_trapz);

S_save = [S_save(B_trapz,:);S_save(C_trapz,:)];
Fnoise_save = [Fnoise_save(B_trapz,:);Fnoise_save(C_trapz,:)];

%% figures setting
Fi = 3;
Fj = Numberof_F/3;

%% smooth tol
p_smooth = [5e-4,5e-4,5e-4,2e-4,...
    5e-4,5e-4,5e-4,5e-4,...
    5e-4,5e-4,5e-4,5e-4];

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
    %% gradients 修改了求导方式！！！！！！！！！！！！
    dsdt = gradient1_s(s,dt10);
    ds2dt2 = gradient2_s(s,dt10);

    lls = [0,find(abs(ds2dt2)>1e-2),discreteT_length];


    %% spline smooth
    F_smooths = csaps(t00,F_noise,p_smooth(i),t00);
     %% filtering and smooth
    X = [];
    dFdt = [];
    dF2dt2 = [];

    F_smooth = [];
    dFdt_smooth = [];
    dF2dt2_smooth = [];
    for k = 1:length(lls)-1
        Xk =  Xk_cut1 + lls(k):lls(k+1)-Xk_cut2;
        X = [X,Xk];
        Length2lls(i,k+1) = Length2lls(i,k) + length(Xk);
%         if isempty(Xk)==0
%             F_smoothk = csaps(t00(Xk),F_noise(Xk),p_smooth(i),t00(Xk));
%             dFdt_smoothk = gradient(F_smoothk,dt10);
%             dF2dt2_smoothk = gradient2(F_smoothk,dt10);
%         end
%         F_smooth = [F_smooth,F_smoothk];
%         dF2dt2 = [dF2dt2, gradient2(F(Xk),dt10)];
%         dF2dt2_smooth = [dF2dt2_smooth, dF2dt2_smoothk];     
%         dFdt = [dFdt, gradient(F(Xk),dt10)];
%         dFdt_smooth = [dFdt_smooth, dFdt_smoothk];
        
        
        F_smoothk = F_smooths(Xk);
        F_smooth = [F_smooth,F_smoothk];
        dFdt_smooth = [dFdt_smooth, gradient1_f(F_smoothk,dt10)];
        dF2dt2_smooth = [dF2dt2_smooth, gradient2_f(F_smoothk,dt10)];     
    end
    
    %% update space
    T_save2 = [T_save2,t10(X)];
    S_save2 = [S_save2,s(X)];

    Fnoise_save2 = [Fnoise_save2,F_noise(X)];
    Fsmooth_save2 = [Fsmooth_save2,F_smooth];
    dFdtsmooth_save2 = [dFdtsmooth_save2, dFdt_smooth];
    dF2dt2smooth_save2 = [dF2dt2smooth_save2,dF2dt2_smooth];

    dsdt_save2 = [dsdt_save2, dsdt(X)];
    ds2dt2_save2 = [ds2dt2_save2, zeros(1,length(X))];

    Length2(i+1) = length(T_save2);
    %% save range
    range_save(i,:) = [min(F_smooth),max(F_smooth),min(s(X)),max(s(X)),...
        min(dFdt_smooth),max(dFdt_smooth),min(dsdt(X)),max(dsdt(X)),...
        min(dF2dt2_smooth),max(dF2dt2_smooth)];
    %% figures
    figure(5);
    subplot(Fi,Fj,i);
    plot(t00,F_noise,'.'); hold on
    plot(t00(X),F_smooth,'.');hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$f$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,s,'--','linewidth',2);hold on
%     ylabel('$s$','Interpreter','latex','Fontsize',18);

    figure(6);
    subplot(Fi,Fj,i);
    plot(t10(X),dFdt_smooth,'.');hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$df/dt$','Interpreter','latex','Fontsize',18);

    figure(7);
    subplot(Fi,Fj,i);
    plot(t10(X),dF2dt2_smooth,'.');hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$d^2f/dt^2$','Interpreter','latex','Fontsize',18); 
end

%%
range_save(Numberof_F+1,:) = [min(range_save(:,1)),max(range_save(:,2)),min(range_save(:,3)),max(range_save(:,4)),...
    min(range_save(:,5)),max(range_save(:,6)),min(range_save(:,7)),max(range_save(:,8)),...
    min(range_save(:,9)),max(range_save(:,10))];

%%
save('exp2-0tp-trapz1.mat','t_scale','Numberof_F','Numberof_Tp','discreteT_length','dt00','t00',...
    'p_smooth','Length2','Length2lls','T_save2',...
    'S_save','S_save2','dsdt_save2','ds2dt2_save2',...
    'Fnoise_save','Fsmooth_save2','dFdtsmooth_save2','dF2dt2smooth_save2',...
    'range_save');