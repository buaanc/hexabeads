function Force_array = contact_law(parameters,delta_array,alphaMax_array)
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
R_from = R;
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
d1 = parameters(5);
d2 = parameters(6);
d3 = parameters(7);
% 	Anorm = pow(alphaNorm,1.137);
e1 = parameters(8);
%  Anorm = (f1*alphaNorm - f2);
f1 = parameters(9);
f2 = parameters(10);
% alphaNorm < g1
g1 = parameters(11);
% Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);
h1 = parameters(12);

%We are returning the vector Force, as long as the vector delta
n_delta = length(delta_array);
Force_array = zeros(n_delta,1);
for i=1:1:n_delta
    
    delta = delta_array(i);
    alphaMax = alphaMax_array(i);
    % alphaMax refers to the level reached in the previous timestep
    % if we are still in the elastic regime, alphaMax = 0
    alphaNorm = alphaMax/alphaY;
    alphaP = c1*alphaMax - c2*alphaY + c3*alphaY*exp(-c4*(alphaNorm-1.0));
    if (alphaP < 0.0)
        alphaP = 0.0;
    end
    if (delta >= alphaP)
        if (delta <= alphaY)
            ke = 4.0/3.0* EStar * sqrt(RStar);

            Force = ke*abs(delta)^1.5;

            %DerivativeForce = 1.5*ke*pow(fabs(delta),0.5);
        else
            if (alphaP == 0)
                % Update the alphaP
                alphaNorm = delta/alphaY;
                alphaP = c1*alphaMax - c2*alphaY + c3*alphaY*exp(-c4*(alphaNorm-1.0));
                if (alphaP<0.0)
                    alphaP = 0.0;
                end


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
                %DerivativeForce = DP*Anorm + DA*pAlpha;

                % Update alphaMax
                alphaMax = delta;
            else
                %Elastic force
                alphaNorm = alphaMax/alphaY;
                % We update alphaP because it changes when we change alphaY
                alphaP = c1*alphaMax - c2*alphaY + c3*alphaY*exp(-c4*(alphaNorm-1.0));
                if (alphaP<0.0)
                    alphaP = 0.0;
                end

                pAlpha = sigmayielding*AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
                if (alphaNorm<g1)
                    Anorm = alphaNorm^e1;
                else
                    Anorm = (f1*alphaNorm - f2);
                end
                FMaxElastic = pAlpha*Anorm;
                Felastic = FMaxElastic*((delta - alphaP)/(alphaMax - alphaP))^h1;

                %Plastic force
                alphaNorm = delta/alphaY;
                alphaP = c1*delta - c2*alphaY + c3*alphaY*exp(-c4*(alphaNorm-1.0));
                if (alphaP<0.0)
                    alphaP = 0.0;
                end

                pAlpha = sigmayielding*AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
                %DP = sigmayielding*AreaY*d2*d3*exp(-d3*(alphaNorm-1.0))/alphaY;
                if (alphaNorm<g1)
                    Anorm = alphaNorm^e1;
                    %DA = 1.0/alphaY*e1*pow(alphaNorm,e1 - 1.0);
                else
                    Anorm = (f1*alphaNorm - f2);
                    %DA = f1/alphaY;
                end
                FPlastic = pAlpha*Anorm;

                %Yield criteria
                if (Felastic < FPlastic)
                    % We need to get back the alphaP corresponding to alphaMax instead of delta
                    alphaNorm = alphaMax/alphaY;
                    alphaP = c1*alphaMax - c2*alphaY + c3*alphaY*exp(-c4*(alphaNorm-1.0));
                    if (alphaP<0.0)
                        alphaP = 0.0;
                    end
                    Force = Felastic;
                    %DerivativeForce = FMaxElastic*h1*1.0/(alphaMax - alphaP)*pow((delta - alphaP)/(alphaMax - alphaP),h1 - 1.0);
                else
                    %Update alphaP
                    alphaNorm = delta/alphaY;
                    alphaP = c1*delta - c2*alphaY + c3*alphaY*exp(-c4*(alphaNorm-1.0));
                    if (alphaP<0.0)
                        alphaP = 0.0;
                    end

                    Force = FPlastic;
                    %DerivativeForce = DP*Anorm + DA*pAlpha;

                    % Update alphaMax
                    alphaMax = delta;
                end
            end
        end
    else
        Force = 0.0;
        %DerivativeForce = 0.0;
    end
    
    Force_array(i) = Force;
end
end