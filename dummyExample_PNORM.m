function dummyExample



function [dke,ke] = stiffness(density_design)
    E_to = 115.00000000000001;
    E_from = 115.00000000000001;
    SIMP_parameter = 4;
    E_to = density_design^SIMP_parameter*E_to;
    nu = 0.3;
    EStar = (E_to*E_from)/(E_to*(1.0-nu^2.0)+E_from*(1-nu^2.0));
    RStar = 4.7629999999999999;
    ke = 4.0/3.0* EStar * sqrt(RStar);
    
    dE = SIMP_parameter*density_design^(SIMP_parameter - 1.0)*E_to;
    dEStar = E_to^2.0*(1-nu^2.0)/(E_from*(1-nu^2.0)+E_to*(1-nu^2.0))^2.0*dE;
    dke = 4.0/3.0* dEStar * sqrt(RStar);
end


    
density_design = 1;
[dk,k] =  stiffness(density_design);

SpaceNorm = 1;
TimeNorm = 1;
m = 3847.2400006626872;
% Number of variables
n = 2;
Tend = 20;
Tspan = [0 Tend]; 

IC = [0 1 0]; % y(t=0) = 1
tol = 1e-5;
initstep = 0.01;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol],'InitialStep',initstep);
[T, Y] = ode23(@(t,y) p_normspring(t,y,k,m,SpaceNorm,TimeNorm),Tspan,IC,options); % Solve ODE

h_todo = diff(T);

fileID = fopen('matlab_timestep.txt','w');
fprintf(fileID,'%12.8f\n',h_todo);
fclose(fileID);

%plot(T, Y(:,2));
%hold on
%plot(T, Y(:,3),'r');
%title('Plot of y as a function of time');
%xlabel('Time'); 
%ylabel('Y(t)');

theta = Y(end,3);
disp('theta')
disp(theta)

timesteps = length(T);
dThetadY = repmat([0;1],1,timesteps);

dFdY = [0, -k/m; 1, 0];
dFdP = [dk*(-Y(:,2))'/m; zeros(1,timesteps)];

A = [0, 0,   0, 0;
    1/2,0,   0, 0;
    0  ,3/4, 0, 0;
    2/9 1/3 4/9 0];
B = [2/9, 1/3, 4/9, 0];
C = [0, 1/2, 3/4, 1];
Stages = 4;

% Finite Difference

dt = 1e-7;

density_design = density_design + dt;
[dk,k] =  stiffness(density_design);
[T_FD, Y_FD] = ode23(@(t,y) p_normspring(t,y,k,m,SpaceNorm,TimeNorm),Tspan,IC,options); % Solve ODE
theta_FD_pre = Y_FD(end,3);


theta_FD = processObj(theta_FD_pre,TimeNorm);
disp('theta_FD')
disp(theta_FD)

density_design = density_design - 2*dt;
[dk,k] =  stiffness(density_design);
[T_FD_2, Y_FD] = ode23(@(t,y) p_normspring(t,y,k,m,SpaceNorm,TimeNorm),Tspan,IC,options); % Solve ODE
theta_FD_2_pre = Y_FD(end,3);

theta_FD_2 = processObj(theta_FD_2_pre,TimeNorm);

[T, T_FD, T_FD_2]
disp('theta_FD_2')
disp(theta_FD_2)

dGdP = (theta_FD - theta_FD_2)/(2*dt);

disp('dGdP')
disp(dGdP)

density_design = density_design + dt;
[dk,k] =  stiffness(density_design);
%Y(:,2)

%length(Y(:,2))


Lambda = zeros(n,timesteps);
mu = zeros(1,timesteps);

TMP_evo = zeros(2,timesteps);
Y2_recalculada = zeros(1,timesteps);
dRdP_S_Total = zeros(timesteps,Stages);
dFdP_S_Total = zeros(timesteps,Stages);
for i=timesteps:-1:1
        %Get primal data
        t = T(i);
        if i>1
            h = - T(i-1) + T(i);
            y = Y(i-1,1:2);
        end
        K_array = zeros(n,Stages);

        Y_i = zeros(n,Stages);
        for istage = 1:Stages  
            Y_i(:,istage) = y';
            Z = zeros(n,1);
            if(istage > 1) 
                for j = 1:(istage-1)
                    Y_i(:,istage) = Y_i(:,istage) + h * A(istage,j)*K_array(:,j);
                end
            end
            %Y_i(:,istage) = Y_i(:,istage) + Z;
            
            T_RK = t + C(istage)*h;
            Y_RK = [Y_i(:,istage); 0];
            K_array_temp = spring(T_RK,Y_RK,k,m,SpaceNorm,TimeNorm);
            K_array(:,istage) = K_array_temp(1:2);           
        end
        
        %Save history
        Y2_recalculada(i) = Y_i(2,1);
        %Evaluate the functions
        
        U_array = zeros(n,Stages);
        V_array = zeros(1,Stages);
        for istage = Stages:-1:1
           dFdY_S =  dFdY;
           %dFdP_S = dFdP(:,i);
           dFdP_S =  [(-dk*Y_i(2,istage))/m; 0];
           dFdP_S_Total(timesteps-i+1,istage) = (-dk*Y_i(2,istage))/m;
           dRdY_S = dthetady_func(t,Y_i(:,istage),k,m,SpaceNorm,TimeNorm);
           dRdP_S =  dthetadP_func(t,Y_i(:,istage),k,m,SpaceNorm,TimeNorm);
           dRdP_S_Total(timesteps-i+1,istage) = dRdP_S;
           TMP = zeros(n,1);
           if istage < Stages              
               for j = (istage + 1):Stages
                   TMP = TMP + h*A(j,istage)*U_array(:,j);
               end     
           end
           TMP = TMP + h*B(istage)*Lambda(:,timesteps-i+1);
           
           U_array(:,istage) = dFdY_S'*TMP + h*B(istage)*dRdY_S;
           V_array(:,istage) = dFdP_S'*TMP + h*B(istage)*dRdP_S;
            
        end
        
        if i>1
            Lambda(:,timesteps-i+2) = Lambda(:,timesteps-i+1) + sum(U_array,2);
            mu(:,timesteps-i+2) = mu(:,timesteps-i+1) + sum(V_array,2);
        end
    
end
disp('Sensitivity')
gradient = mu(end)^(1 / TimeNorm);
disp(gradient)

%mu = mu' + 1/(Tend)*T;
figure
plot(T,mu)
hold on
plot(T,Lambda,'--')
end



function dy = spring(t,y,k,m,SpaceNorm,TimeNorm)

dy = zeros(3,1);    % a column vector

dy(1) = 1/m*(-k*y(2));
dy(2) = y(1);

dy(3) = y(2);


end

function Obj = processObj(Objpre,TimeNorm)

    Obj = Objpre^(1 / TimeNorm);

end

function dy = p_normspring(t,y,k,m,SpaceNorm,TimeNorm)

dy = zeros(3,1);    % a column vector

dy(1) = 1/m*(-k*y(2));
dy(2) = y(1);

Force = k*y(2); 

Force_obj = (Force^SpaceNorm)^(TimeNorm / SpaceNorm);
dy(3) = Force_obj;


end

function dthetady = dthetady_func(t,y,k,m,SpaceNorm,TimeNorm)

Deriv_Force = k; 
dthetady = zeros(2,1);

Force = k*y(2); 

deriv_coeff = TimeNorm / SpaceNorm*(Force^SpaceNorm)^(TimeNorm / SpaceNorm - 1 )*...
    (SpaceNorm * Force^(SpaceNorm - 1));
dthetady(2)= deriv_coeff*Deriv_Force;


end

function dthetadP = dthetadP_func(t,y,k,m,SpaceNorm,TimeNorm)
density_design = 1;
E_to = 115.00000000000001;
E_from = 115.00000000000001;
SIMP_parameter = 4;
E_to = density_design^SIMP_parameter*E_to;
nu = 0.3;
EStar = (E_to*E_from)/(E_to*(1.0-nu^2.0)+E_from*(1-nu^2.0));
RStar = 4.7629999999999999;
ke = 4.0/3.0* EStar * sqrt(RStar);

dE = SIMP_parameter*density_design^(SIMP_parameter - 1.0)*E_to;
dEStar = E_to^2.0*(1-nu^2.0)/(E_from*(1-nu^2.0)+E_to*(1-nu^2.0))^2.0*dE;
dke = 4.0/3.0* dEStar * sqrt(RStar);
    
Deriv_Force = dke*y(2); 

Force = k*y(2); 

deriv_coeff = TimeNorm / SpaceNorm*(Force^SpaceNorm)^(TimeNorm / SpaceNorm - 1 )*...
    (SpaceNorm * Force^(SpaceNorm - 1));
dthetadP= deriv_coeff*Deriv_Force;


end