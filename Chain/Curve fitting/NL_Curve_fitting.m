% --------------------------- Read experimental data --------------------

%Grabs the smoothed experimental data and calculates the state variable
%alphaMax, call it with 1 (reading_data(1)) to understand what is doing,
%the new_delta, new_force and new_alphaMax are just the same data but
%cropped. I discard the data of the elastic reloading.
[delta, alphaMax, force,...
    new_delta,new_force, new_alphaMax,...
    alphaMaxUnloading, delta_unloaded,...
    delta_plastic_regime, force_plastic_regime, ...
    delta_elastic_regime, force_elastic_regime, alphaMax_elastic_regime,...
    smoothed_delta, smoothed_force]...
    = reading_data(0);


% ------------- Fitting parameters plastic regime --------------
% pAlpha = sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
d1 = 2.48;
Initial_parameters_plastic(1) = d1;
d2 = 1.41;
Initial_parameters_plastic(2) = d2;
d3 = 0.098;
Initial_parameters_plastic(3) = d3;
% 	Anorm = pow(alphaNorm,e1);
e1 =1.137;
Initial_parameters_plastic(4) = e1;
%  Anorm = (f1*alphaNorm - f2);
f1 = 2.37076627;
Initial_parameters_plastic(5) = f1;
f2 = 59.955374401;
Initial_parameters_plastic(6) = f2;
% alphaNorm < g1
g1 = 177.57;
Initial_parameters_plastic(7) = g1;

fitting_func_plastic_regime = @(parameters,delta_plastic_regime)...
                        contact_law_plastic_regime(parameters,delta_plastic_regime);
opts = optimset('MaxFunEvals',1e7,'MaxIter',1e6,'Display','Iter','PlotFcns',@optimplotx);
x = lsqcurvefit(fitting_func_plastic_regime,Initial_parameters_plastic,...
                delta_plastic_regime,force_plastic_regime,[],[],opts);
            
fitted_param_plastic = x;

% Checking the fitting

fitting_force_plastic_regime = contact_law_plastic_regime(x,delta_plastic_regime);
figure
plot(delta_plastic_regime,fitting_force_plastic_regime);
hold on
plot(delta_plastic_regime,force_plastic_regime,'r-');
xlabel('Strain')
ylabel('Force')
title('Comparing plastic regime')

% Now we have left only the unloading/reloading curves
% Convert the cell arrays into a continuous array


delta_elastic = cell2mat(delta_elastic_regime);
force_elastic = cell2mat(force_elastic_regime);
alphaMax_elastic = cell2mat(alphaMax_elastic_regime);

% Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);
Initial_parameters_elastic =zeros(5,1);

c1 = 0.9981927511;
Initial_parameters_elastic(1) = c1;
%c2 = 26.934752932; 
c2 = 25.934752932;
Initial_parameters_elastic(2) = c2;
c3 = 24.9865601805;
Initial_parameters_elastic(3) = c3;
%c4 = 0.0015; 
c4 = 0.015;
Initial_parameters_elastic(4) = c4;
h1 = 1.35;

Initial_parameters_elastic(5) = h1;

fitting_func_elastic_regime = @(parameters,delta_elastic)...
                        contact_law_elastic_regime(parameters,delta_elastic,alphaMax_elastic,...,
                        fitted_param_residual, fitted_param_plastic);

opts = optimset('MaxFunEvals',1e7,'MaxIter',1e6,'Display','Iter','PlotFcns',@optimplotx);

x = lsqcurvefit(fitting_func_elastic_regime,Initial_parameters_elastic,...
                delta_elastic,force_elastic,[],[],opts);

            
            
% Checking the fitting

fitting_force_elastic_regime = contact_law_elastic_regime(x,delta_elastic,alphaMax_elastic,...,
                        fitted_param_residual, fitted_param_plastic);
figure
plot(delta_elastic,fitting_force_elastic_regime);
hold on
plot(delta_elastic,force_elastic,'r-');
xlabel('Strain')
ylabel('Force')
title('Comparing elastic regime')            
          
            

% ------------ Plotting to check the overall fitting --------------------
% fitted_param are the ones that fit the experimental values,
% raj_fitted_param are the original ones, raj_bool pick either fitted_param
% (0) or raj_fitted_param (1)

all_parameters = [x(1:4); fitted_param_plastic'; x(5)];
Force_Raj_fit = contact_law(all_parameters,delta,alphaMax);
                
figure          
plot(delta,Force_Raj_fit);
hold on
plot(smoothed_delta,smoothed_force,'r-');
xlabel('Strain')
ylabel('Force')
title('Comparing contact laws after fitting')


                
