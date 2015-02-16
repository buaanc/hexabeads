function ClipData

SS=csvread('Yall.csv',1,0); % Case A: Loading-unloading data from 1 to 6 kN in 1 kN increments
DD=csvread('Dall.csv',1,0); % Case B: Loading-unloading data from 1.5 to 5.5 kN in 1 kN increments
AA=csvread('All2500.csv',1,0); % Case C: Load to and unload from 2.5 kN
BB=csvread('All5kN.csv',1,0); % Case D: Load to and unload from 5 kN
CC=csvread('All10kN.csv',1,0); % Case E: Load to and unload from 10 kN


Load=SS(:,2);Strain=SS(:,3); % Load and Strain from Case A
LD=DD(:,2);SD=DD(:,3); % Load and Strain from Case B
LA=AA(:,1);SA=AA(:,2); %"" Case C
LB=BB(:,1);SB=BB(:,2); %"" Case D
LC=CC(:,1);SC=CC(:,2); % "" Case E

% Offsets
S1=-.041; % Offset for Case E
S2=-.047; % offset for Case D
S3=-.0366; % offset for Case B

figure
hold on
plot(-1.*Strain,-1.*Load,'k');
plot(-1.*(SD+S3),-1.*LD,'r+')
plot(-1.*SA,-1.*LA,'b')
plot(-1.*(SB+S2),-1.*LB,'m')
plot(-1.*(SC+S1),-1.*LC,'g')

hold off
xlabel('Deformation (mm)','FontSize',16)
ylabel('Load (kN)','FontSize',16)
legend('Case A','Case B','Case C','Case D','Case E')
grid on
set(gca,'FontSize',14)

% I'm going to split this case and get rid of the second part
case_B_strain = -1.*(SD+S3);
case_B_force = -1.*LD;

% Look up the gap where there's no unloading curve. It is where there's a
% gap in the data.
difference = diff(case_B_force);
[maximo,new_start] = max(abs(difference));

case_B_1_strain = case_B_strain(1:new_start);
case_B_1_force = case_B_force(1:new_start);

case_B_2_strain = case_B_strain(new_start+1:end);
case_B_2_force = case_B_force(new_start+1:end);


Experimental_data = cell(5,1);
Experimental_data{1} = [-1.*Strain, -1.*Load];
Experimental_data{2} = [case_B_1_strain, case_B_1_force];
%Experimental_data{3} = [case_B_2_strain, case_B_2_force];
Experimental_data{3} = [-1.*SA, -1.*LA];
Experimental_data{4} = [-1.*(SB+S2), -1.*LB];
Experimental_data{5} = [-1.*(SC+S1), -1.*LC];

% Plotting to 1 to see the whole process to get the meaningful data
plotting = 0;

load('parameters.mat')
 
for j=1:5
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
