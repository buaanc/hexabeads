function Force_array = contact_law_plastic_regime(parameters,delta_array)
% ----------------- Material parameters -----------------------
% Young Modulus
E = 115; %Gpa
nu = 0.3;
% Same material for both beads
E_to = E;
E_from = E;
% Combined young modulus
EStar = (E_to*E_from)/(E_to*(1.0-nu^2)+E_from*(1-nu^2));

% Radius
R = 4.763; % mm
% Same radius
R_to = R;
R_from = 1.79;
%Combined radius
RStar = (R_to*R_from)/(R_to+R_from);

%Yield stress
sigmayielding = 0.55; %Gpa

% Yield parameters
alphaY = (pi/2)^2* (sigmayielding*1.6/EStar)^2 * RStar;
AreaY = pi*RStar*alphaY;

% ------------------------ Fitting parameters --------------------------



% pAlpha = sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
d1 = parameters(1);
d2 = parameters(2);
d3 = parameters(3);
% 	Anorm = pow(alphaNorm,1.137);
e1 = parameters(4);
%  Anorm = (f1*alphaNorm - f2);
f1 = parameters(5);
f2 = parameters(6);
% alphaNorm < g1
g1 = parameters(7);


tol = 1e-10;

%We are returning the vector Force, as long as the vector delta
n_delta = length(delta_array);
Force_array = zeros(n_delta,1);

for i=1:1:n_delta
    
    delta = delta_array(i);
    alphaNorm = delta/alphaY;
    pAlpha =sigmayielding*AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
    %DP = sigmayielding*AreaY*d2*d3*exp(-d3*(alphaNorm-1.0))/alphaY;
    if (alphaNorm<g1)
        Anorm = alphaNorm^e1;
        %DA = 1.0/alphaY*e1*pow(alphaNorm,e1 - 1.0);
    else
        Anorm = (f1*alphaNorm - f2);
        %DA = f1/alphaY;
    end
    Force = pAlpha*Anorm;
    Force_array(i) = Force;
end