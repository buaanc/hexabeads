close all; 
clear all;

CaseA=csvread('Yall.csv',1,0); % Case A: Loading-unloading data from 1 to 6 kN in 1 kN increments
CaseB=csvread('Dall.csv',1,0); % Case B: Loading-unloading data from 1.5 to 5.5 kN in 1 kN increments
CaseC=csvread('All2500.csv',1,0); % Case C: Load to and unload from 2.5 kN
CaseD=csvread('All5kN.csv',1,0); % Case D: Load to and unload from 5 kN
CaseE=csvread('All10kN.csv',1,0); % Case E: Load to and unload from 10 kN


Load_CaseA=CaseA(:,2);Deformation_CaseA=CaseA(:,3); % Load and Deformation from Case A
Load_CaseB=CaseB(:,2);Deformation_CaseB=CaseB(:,3); % Load and Deformation from Case B
Load_CaseC=CaseC(:,1);Deformation_CaseC=CaseC(:,2); %"" Case C
Load_CaseD=CaseD(:,1);Deformation_CaseD=CaseD(:,2); %"" Case D
Load_CaseE=CaseE(:,1);Deformation_CaseE=CaseE(:,2); % "" Case E



figure
hold on
plot(Deformation_CaseA,Load_CaseA,'k');
plot(Deformation_CaseB,Load_CaseB,'r')
plot(Deformation_CaseC,Load_CaseC,'b')
plot(Deformation_CaseD,Load_CaseD,'m')
plot(Deformation_CaseE,Load_CaseE,'g')

hold off
xlabel('Deformation (mm)','FontSize',16)
ylabel('Load (kN)','FontSize',16)
grid on
legend('A','B','C','D','E')
set(gca,'FontSize',14)


xlim([0 .5])
ylim([0 10.1])
