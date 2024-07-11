clc; clear; close all;
iter = 100;
S = [25,50,75,100,125,150,175,200,225,250];
load('nominalFDR_005.txt')
load('nominalFDR_01.txt')
load('FDR_005.txt')
load('FDR_01.txt')
load('stdFDR_005.txt')
load('stdFDR_01.txt')
load('POWER_005.txt')
load('POWER_01.txt')
load('stdPOWER_005.txt')
load('stdPOWER_01.txt')
%% plot errorbar
subplot(3,1,1)
errorbar(S,POWER_005,stdPOWER_005./sqrt(iter),'b-o','LineWidth',1)
hold on
errorbar(S,POWER_01,stdPOWER_01./sqrt(iter),'r-s','LineWidth',1)
hold on
legend('q=0.05','q=0.1','','')
grid on
axis([20,255,0.97,1.01]);
ylabel('TP')
subplot(3,1,2:3)
errorbar(S,FDR_005,stdFDR_005./sqrt(iter),'b-o','LineWidth',1)
hold on
errorbar(S,FDR_01,stdFDR_01./sqrt(iter),'r-s','LineWidth',1)
hold on
plot(S,nominalFDR_005,'black--',S,nominalFDR_01,'black--','LineWidth',1)
hold on
legend('q=0.05','q=0.1','','')
grid on
axis([20,255,0.03,0.1]);
xlabel('Sparsity');
ylabel('TgFDR')