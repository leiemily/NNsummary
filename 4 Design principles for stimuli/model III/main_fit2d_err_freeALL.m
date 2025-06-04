clear;
%% main
j = 1;
is1 = 9;
is2 = 10;
% it1 = 12;
% it2 = 16;
it1 = 12;
it2 = 10;
il1 = 15;
il2 = 16;

err_mean = zeros(1,5);

%% paths
folderPath_NNdata = '.\dataset2\';

%% L
setName1 = ['bias2-',num2str(j),'tp-test2d-freeL'];
load(fullfile(folderPath_NNdata, setName1),'err1','err2','err3');
err_mean(1) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% LS
setName2 = ['bias2-',num2str(j),'tp-test2d-freeLS'];
load(fullfile(folderPath_NNdata, setName2),'err1','err2','err3');
err_mean(2) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% LT
setName3 = ['bias2-',num2str(j),'tp-test2d-freeLT'];
load(fullfile(folderPath_NNdata, setName3),'err1','err2','err3');
err_mean(3) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% TSL
setName4 = ['bias2-',num2str(j),'tp-test2d-freeTSL'];
load(fullfile(folderPath_NNdata, setName4),'err1','err2','err3');
err_mean(4) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;

%% TS
setName5 = ['bias2-',num2str(j),'tp-test2d-freeTS'];
load(fullfile(folderPath_NNdata, setName5),'err1','err2','err3');
err_mean(5) = (err1(end,is1) + err1(end,is2) + err2(end,it1) + err2(end,it2) + err3(end,il1) + err3(end,il2))/6;