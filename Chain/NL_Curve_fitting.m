% --------------------------- Read experimental data --------------------

%Grabs the smoothed experimental data and calculates the state variable
%alphaMax, call it with 1 (reading_data(1)) to understand what is doing,
%the new_delta, new_force and new_alphaMax are just the same data but
%cropped. I discard the data of the elastic reloading.
[delta, alphaMax, force,new_delta,new_force, new_alphaMax] = reading_data(0);

% ------------------------ Fitting parameters --------------------------
Initial_parameters = zeros(1,12);
% 	alphaP = c1*delta - c2*alphaY + c3*alphaY*exp(-c4*(delta/alphaY-1.0));
c1 = 1.022; 
c1 = 0.9981927511;
Initial_parameters(1) = c1;
c2 = 26.934752932; 
c2 = 25.934752932;
Initial_parameters(2) = c2;
c3 = 24.9865601805;
Initial_parameters(3) = c3;
c4 = 0.0015; 
c4 = 0.015;
Initial_parameters(4) = c4;
% pAlpha = sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
d1 = 3.45;
d1 = 2.48;
Initial_parameters(5) = d1;
d2 = d1 - 1.07; %1.41;
Initial_parameters(6) = d2;
d3 = 0.098;
Initial_parameters(7) = d3;
% 	Anorm = pow(alphaNorm,1.137);
e1 = 1.037; 
e1 =1.137;
Initial_parameters(8) = e1;
%  Anorm = (f1*alphaNorm - f2);
f1 = 2.77076627; 
f1 = 2.37076627;
Initial_parameters(9) = f1;
f2 = 59.955374401;
Initial_parameters(10) = f2;
% alphaNorm < g1
g1 = 243.57;
g1 = 177.57;
Initial_parameters(11) = g1;
% Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);
h1 = 1.35;
Initial_parameters(12) = 1.35;

% ------------ Plotting to check the initial fitting --------------------
Force_Raj = contact_law(Initial_parameters,delta,alphaMax);

plot(delta,Force_Raj);
hold on
plot(delta,force,'r-');
hold on
plot(new_delta,new_force,'k')
xlabel('Strain')
ylabel('Force')
title('Comparing contact laws')

%---------------------Non Linear Curve fitting ------------------------
fitting_func = @(parameters,new_delta)contact_law(parameters,new_delta,new_alphaMax);
opts = optimset('MaxFunEvals',1e7,'MaxIter',1e6,'Display','Iter','PlotFcns',@optimplotx);
x = lsqcurvefit(fitting_func,Initial_parameters,new_delta,new_force,[],[],opts);

% ------------ Plotting to check the initial fitting --------------------
Force_Raj_fit = contact_law(x,delta,alphaMax);

plot(delta,Force_Raj_fit);
hold on
plot(delta,force,'r-');
xlabel('Strain')
ylabel('Force')
title('Comparing contact laws after fitting')