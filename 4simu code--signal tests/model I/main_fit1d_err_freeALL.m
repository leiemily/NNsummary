clear;
%% main
j = 1;
% is1 = 15;
% is2 = 16;
% it1 = 15;
% it2 = 16;
is1 = 11;
is2 = 12;
it1 = 10;
it2 = 13;
il1 = 15;
il2 = 16;

err_mean = zeros(1,6);

%% paths
folderPath_NNdata = '.\dataset2\';

%% T
setName1 = ['frac1-',num2str(j),'tp-test1d-freeT'];
load(fullfile(folderPath_NNdata, setName1),'err1','err2','err3');
err_mean(1) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% L
setName2 = ['frac1-',num2str(j),'tp-test1d-freeL'];
load(fullfile(folderPath_NNdata, setName2),'err1','err2','err3');
err_mean(2) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% TS
setName3 = ['frac1-',num2str(j),'tp-test1d-freeTS'];
load(fullfile(folderPath_NNdata, setName3),'err1','err2','err3');
err_mean(3) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% LS
setName4 = ['frac1-',num2str(j),'tp-test1d-freeLS'];
load(fullfile(folderPath_NNdata, setName4),'err1','err2','err3');
err_mean(4) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% LT
setName5 = ['frac1-',num2str(j),'tp-test1d-freeLT'];
load(fullfile(folderPath_NNdata, setName5),'err1','err2','err3');
err_mean(5) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% TSL
setName6 = ['frac1-',num2str(j),'tp-test1d-freeTSL'];
load(fullfile(folderPath_NNdata, setName6),'err1','err2','err3');
err_mean(6) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;