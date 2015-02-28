function Force_array = contact_law_elastic_regime(parameters,delta_array,alphaMax_array,fitted_param_plastic)
% ----------------- Material parameters -----------------------
% Young Modulus
E = 115; %Gpa
nu = 0.3;
% Same material for both beads
E_to = E;
E_from = 200;
% Combined young modulus
EStar = (E_to*E_from)/(E_to*(1.0-nu^2)+E_from*(1-nu^2));

% Radius
R = 4.763; % mm
% Same radius
R_to = R;
R_from = (sqrt(2.0)-1.0)*R;
%Combined radius
RStar = (R_to*R_from)/(R_to+R_from);


%Yield stress
sigmayielding = 0.55; %Gpa

% Yield parameters
alphaY = (pi/2)^2* (sigmayielding*1.6/EStar)^2 * RStar;
AreaY = pi*RStar*alphaY;

% ------------------------ Fitting parameters --------------------------

% 	alphaP = c1*delta - c2*alphaY + c3*alphaY*exp(-c4*(delta/alphaY-1.0));
c1 = parameters(1);
c2 = parameters(2);
c3 = parameters(3);
c4 = parameters(4);

% pAlpha = sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
d1 = fitted_param_plastic(1);
d2 = fitted_param_plastic(2);
d3 = fitted_param_plastic(3);
% 	Anorm = pow(alphaNorm,1.137);
e1 = fitted_param_plastic(4);
%  Anorm = (f1*alphaNorm - f2);
f1 = fitted_param_plastic(5);
f2 = fitted_param_plastic(6);
% alphaNorm < g1
g1 = fitted_param_plastic(7);

% Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);
h1 = parameters(5);

%We are returning the vector Force, as long as the vector delta
n_delta = length(delta_array);
Force_array = zeros(n_delta,1);

%alphaP = c1*alphaMax_array - c2*alphaY + c3*alphaY*exp(-c4*(alphaMax_array/alphaY-1.0));

for i=1:1:n_delta

    delta = delta_array(i);
    alphaMax = alphaMax_array(i);
    %Elastic force
    alphaNorm = alphaMax/alphaY;
    % We update alphaP because it changes when we change alphaY
    alphaP = c1*alphaMax - c2*alphaY + c3*alphaY*exp(-c4*(alphaNorm-1.0));
    if (alphaP<0.0) 
        alphaP = 0.0;
    end
    if (alphaP > delta)
        alphaP = delta;
    end
    
    if alphaP > alphaMax
        a= 1;
    end

    pAlpha = sigmayielding*AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
    if (alphaNorm<g1)
        Anorm = alphaNorm^e1;
    else
        Anorm = (f1*alphaNorm - f2);
    end
    FMaxElastic = pAlpha*Anorm;
    Force = FMaxElastic*((delta - alphaP)/(alphaMax - alphaP))^h1;

    Force_array(i) = Force;
end

end