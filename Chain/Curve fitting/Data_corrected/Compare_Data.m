function ClipData


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
legend('Case A','Case B','Case C','Case D','Case E')
grid on
set(gca,'FontSize',14)

Experimental_data = cell(5,1);
Experimental_data{1} = [Deformation_CaseA, Load_CaseA];
Experimental_data{2} = [Deformation_CaseB, Load_CaseB];
Experimental_data{3} = [Deformation_CaseC, Load_CaseC];
Experimental_data{4} = [Deformation_CaseD, Load_CaseD];
Experimental_data{5} = [Deformation_CaseE, Load_CaseE];

% Plotting to 1 to see the whole process to get the meaningful data
plotting = 1;

load('parameters.mat')
 
figure
for j=1:1
    [delta, alphaMax, force,...
        new_delta,new_force, new_alphaMax, ...
        alphaMaxUnloading,delta_unloaded,...
        delta_plastic_regime, force_plastic_regime,...
        delta_elastic_regime, force_elastic_regime, alphaMax_elastic_regime,...
        smoothed_delta, smoothed_force] = validate_data(plotting,Experimental_data{j},j);



    Force_Raj_fit = contact_law(all_parameters,delta,alphaMax);

    figure 
    plot(delta,Force_Raj_fit);
    hold on
    plot(smoothed_delta,smoothed_force,'r-');
    xlabel('Strain')
    ylabel('Force')
    title('Comparing contact laws after fitting')
end

end
