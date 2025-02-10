clear;
%% mixsin1 先全局光滑再分段求导 200s
load('tumbling_fraction_mixsinall_dt05');
Snoise_save = S_save;
num_fs = [5,6,7];
Snoise_save = Snoise_save(num_fs,:);
Fnoise_save = Fnoise_save(num_fs,:);

Numberof_F = size(Fnoise_save,1);

s_cut = 4:length(t)-1;
Fnoise_save = Fnoise_save(:,s_cut);
Snoise_save = Snoise_save(:,s_cut);

%% rescaling
t_scale = 100;
dt00 = 0.5;
t00 = t(s_cut)+ 0.1 -2;%改成0:...
discreteT_length = length(t00);

%% para
pf_smooth = 5e-4*ones(1,Numberof_F);
ps_smooth = 5e-4*ones(1,Numberof_F);
%% noise related
S_save = zeros(Numberof_F,discreteT_length);
Fsmooth_save = zeros(Numberof_F,discreteT_length);

%% iter
for i = 1:Numberof_F
    s_noise = Snoise_save(i,:);
    F_noise = Fnoise_save(i,:);
    %% smooth
    F_smooth = csaps(t00,F_noise,pf_smooth(i),t00);
    S_smooth = csaps(t00,s_noise,ps_smooth(i),t00);
    %% update space
    S_save(i,:) = S_smooth;
    Fsmooth_save(i,:) = F_smooth;
end

% Fi = Numberof_F;
% Fj = 1;
% for i = 1:Numberof_F
%     figure(13);   
%     subplot(Fi,Fj,i);
%     plot(t00,Fnoise_save(i,:),'.'); hold on
%     plot(t00,Fsmooth_save(i,:),'-','linewidth',2); hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$f$','Interpreter','latex','Fontsize',18); 
%     yyaxis right;
%     plot(t00,Snoise_save(i,:),'--','linewidth',2);hold on
%     plot(t00,S_save(i,:),'-','linewidth',2);hold on
%     ylabel('$s$','Interpreter','latex','Fontsize',18); 
% 
% end

%% 分段200s
num_tpart = 6;
window_length = 8;
part_length = 406;
s_save = [];
snoise_save = [];
fnoise_save = [];
fsmooth_save = [];
for num = 1:num_tpart
    t_part = part_length*(num-1)-(window_length*(num-1)-1):part_length*num -window_length*(num-1);
    s_save = [s_save;S_save(:,t_part)];
    snoise_save = [snoise_save;Snoise_save(:,t_part)];
    fnoise_save = [fnoise_save;Fnoise_save(:,t_part)];
    fsmooth_save = [fsmooth_save;Fsmooth_save(:,t_part)];
end

%% 重排，来自同一条的放一起
num_rearrange = [];
for i = 1:length(num_fs)
    for num = 1:num_tpart
        num_rearrange = [num_rearrange,i+(num-1)*length(num_fs)];
    end
end

S_save = s_save(num_rearrange,:);
Snoise_save = snoise_save(num_rearrange,:);
Fnoise_save = fnoise_save(num_rearrange,:);
Fsmooth_save = fsmooth_save(num_rearrange,:);
Numberof_F = size(Fsmooth_save,1);

%% 调整顺序
% C_mix = [17,18,21,22,...
%     2,4,5,6,7,8,9,10];
C_mix = [9,2,3,11,4,5];
A_mix = 1:Numberof_F;
B_mix = setdiff(A_mix,C_mix);

S_save = [S_save(B_mix,:);S_save(C_mix,:)];
Snoise_save = [Snoise_save(B_mix,:);Snoise_save(C_mix,:)];
Fnoise_save = [Fnoise_save(B_mix,:);Fnoise_save(C_mix,:)];
Fsmooth_save = [Fsmooth_save(B_mix,:);Fsmooth_save(C_mix,:)];
%% time scaling
t00 = 0:dt00:(part_length-1)*dt00;
dt10 = dt00/t_scale;
discreteT_length = length(t00);

%% 
dsdt_save = zeros(Numberof_F,discreteT_length);
ds2dt2_save = zeros(Numberof_F,discreteT_length);
dFdtsmooth_save = zeros(Numberof_F,discreteT_length);
dF2dt2smooth_save = zeros(Numberof_F,discreteT_length);

%% iter
for i = 1:Numberof_F
    s = S_save(i,:);
    F_smooth = Fsmooth_save(i,:);
    %% gradients
    dsdt = gradient1_s(s,dt10);
    dFdt_smooth = gradient1_f(F_smooth,dt10);
    
    ds2dt2 = gradient2_s(s,dt10);
    dF2dt2_smooth = gradient2_f(F_smooth,dt10);  
    %% update space
    dsdt_save(i,:) = dsdt;
    ds2dt2_save(i,:) = ds2dt2;
    
    dFdtsmooth_save(i,:) = dFdt_smooth;
    dF2dt2smooth_save(i,:) = dF2dt2_smooth;
    
end
%% s_smooth minus, cut
s_cut = 5:length(t00)-1;
S_save = S_save(:,s_cut);
Snoise_save = Snoise_save(:,s_cut);
Fsmooth_save = Fsmooth_save(:,s_cut);
Fnoise_save = Fnoise_save(:,s_cut);
dsdt_save = dsdt_save(:,s_cut);
ds2dt2_save = ds2dt2_save(:,s_cut);
dFdtsmooth_save = dFdtsmooth_save(:,s_cut);
dF2dt2smooth_save = dF2dt2smooth_save(:,s_cut);

t00 = t00(s_cut)-t00(5);
t10 = t00/t_scale;
discreteT_length = length(t00);
%% save range
range_save = [min(Fsmooth_save,[],2),max(Fsmooth_save,[],2),min(S_save,[],2),max(S_save,[],2),...
    min(dFdtsmooth_save,[],2),max(dFdtsmooth_save,[],2),min(dsdt_save,[],2),max(dsdt_save,[],2),...
    min(dF2dt2smooth_save,[],2),max(dF2dt2smooth_save,[],2),min(ds2dt2_save,[],2),max(ds2dt2_save,[],2)];
range_save(Numberof_F+1,:) = [min(range_save(:,1)),max(range_save(:,2)),min(range_save(:,3)),max(range_save(:,4)),...
    min(range_save(:,5)),max(range_save(:,6)),min(range_save(:,7)),max(range_save(:,8)),...
    min(range_save(:,9)),max(range_save(:,10)),min(range_save(:,11)),max(range_save(:,12))];

%% % figures
Fi = length(num_fs);
Fj = Numberof_F/Fi;
for i = 1:Numberof_F
    figure(13);   
    subplot(Fi,Fj,i);
    plot(t00,Fnoise_save(i,:),'.'); hold on
    plot(t00,Fsmooth_save(i,:),'-','linewidth',2); hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$f$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t00,Snoise_save(i,:),'--');hold on
    plot(t00,S_save(i,:),'-','linewidth',2);hold on
%     ylabel('$s$','Interpreter','latex','Fontsize',18); 

    figure(4);   
    subplot(Fi,Fj,i);
    plot(t10,dFdtsmooth_save(i,:),'--','linewidth',2);hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$df/dt$','Interpreter','latex','Fontsize',18);
    yyaxis right;
    plot(t10,dsdt_save(i,:),':','linewidth',2);
%     ylabel('$ds/dt$','Interpreter','latex','Fontsize',18); 

    figure(5);   
    subplot(Fi,Fj,i);
    plot(t10,dF2dt2smooth_save(i,:),'--','linewidth',2);hold on
%     xlabel('$t$','Interpreter','latex','Fontsize',18);
%     ylabel('$d^2f/dt^2$','Interpreter','latex','Fontsize',18); 
    yyaxis right;
    plot(t10,ds2dt2_save(i,:),':','linewidth',2);
%     ylabel('$d^2s/dt^2$','Interpreter','latex','Fontsize',18); 
end
%%
save('exp2-0tp-mix1.mat','t_scale','Numberof_F','discreteT_length','dt00','t00',...
    'pf_smooth','ps_smooth',...
    'S_save','dsdt_save','ds2dt2_save',...
    'Fnoise_save','Fsmooth_save','dFdtsmooth_save','dF2dt2smooth_save',...
    'range_save');