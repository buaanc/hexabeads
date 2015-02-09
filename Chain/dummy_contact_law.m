function dummy

k1 = 1;
k2 = 2;
h = 0.1;
scale_factor = 10000;

format long

t_total = 6.206;
%t_total = 4.206;
t_total = 1;
[obj_function_1,obj_function_ultimo,delta_history,alphaMaxhistory,delta_intermedio ] = primal_problem(k1,k2,t_total,h,scale_factor);


obj_function_1

direct_diff(k1,k2,t_total,h,scale_factor)
% 
% 
sensitivity = ...
      adjoint_method (k1,k2,t_total,h, delta_history, alphaMaxhistory,scale_factor,delta_intermedio);
% 
% disp('Function analitica')
% phi = 2 / k1^2 * cos(sqrt(k1) * t_total) + t_total^2 /k1 - 2 / k1^2;
% disp(phi)
% 
% %t_analytic = 0:h:t_total;
% %solution = -2 / (k1 * sqrt(k1)) * sin(sqrt(k1) * t) +
% 
% disp('Sensitivity analitica')
% dphi = -4 / k1^3 * cos(sqrt(k1) * t_total) -2/k1^2*sin(sqrt(k1)*t_total) * t_total/(2*sqrt(k1)) + 4 / k1^3 ...
%     - t_total^2/k1^2;
% disp(dphi)
% 
% 
% disp('Sensitivity analitica ultimo')
% dphi = 3/ k1^(5/2) * sin(sqrt(k1) * t_total) -2/k1^(3/2)*cos(sqrt(k1)*t_total) * t_total/(2*sqrt(k1)) - 2*t_total/k1^2;
% disp(dphi)

%Finite difference
dt = 1e-8;
k1 = k1 + dt;
[obj_function_2,obj_function_ultimo_2,delta_history, alphaMaxhistory,delta_intermedio]  = primal_problem(k1,k2,t_total,h,scale_factor);


k1 = k1 - 2*dt;
[obj_function_3, obj_function_ultimo_3, delta_history,alphaMaxhistory,delta_intermedio] = primal_problem(k1,k2,t_total,h,scale_factor);
k1 = k1 + dt;

sensitivity_FD = (obj_function_2 - obj_function_3)/(2*dt)

%sensitivity_FD_ultimo = (obj_function_ultimo_2 - obj_function_ultimo_3)/(2*dt)




function [obj_function, obj_function_ultimo, delta, alphaMaxhistory,delta_intermedio] = ...
                    primal_problem(k1,k2,t_total,h,scale_factor)
    
    alphaMax = 0;
    t = 0;
    obj_function = 0;

    m = 1;



    y = zeros(2,1);

    N_steps = length(0:h:t_total);
    delta = zeros(1,N_steps);
    delta_intermedio = zeros(1,N_steps);
    alphaMaxhistory = zeros(1,N_steps);
    
    
    t_array = zeros(1,N_steps);
    force_array = zeros(1,N_steps);
    k = 1;
    plast = 0;
    for i=0:h:t_total-h
        delta(k) = y(2);
        %Update alphaMax with previous value
        if plast == 1
            alphaMax = y(2);
        end    
        alphaMaxhistory(k) = alphaMax;
        % Midpoint method
        
        % First evaluation
        [force, plast,int_force] = RHS(y(2),alphaMax,k1,k2,m,t);

        %obj_function = obj_function + 1/2*h*y(2)/scale_factor;
        y(2) = y(2) + 1/2*h*y(1);
        y(1) = y(1) + 1/2*h*force;
        
        % Save stage variable
        delta_intermedio(k) = y(2);
        
        
        % Value we'll feed to the next stage
        t = t + 1/2*h; 
        
        % Call again the RHS
        [force, plast,int_force] = RHS(y(2),alphaMax,k1,k2,m,t);
        
        obj_function = obj_function + h*y(2)/scale_factor;
        y(2) = y(2) + h*y(1);
        y(1) = y(1) + h*force;
        
        t = t + 1/2*h; 
        

%         % Call again the RHS
%         [force, plast,int_force] = RHS(y(2),alphaMax,k1,k2,m,t);
        
        t_array(k) = t;
        force_array(k) = int_force;
        
    
        
        k = k + 1;

        
    end

%     plot(delta,force_array)
%     xlabel('delta')
%     ylabel('force')
%     
%     figure(2)
%     plot(t_array,alphaMaxhistory)
    
    
    obj_function_ultimo = y(2);
    
    %Save last values of delta and alphaMax
    delta(k) = y(2);
    if plast == 1
            alphaMax = y(2);
    end 
    alphaMaxhistory(k) = alphaMax;
    
%     figure
%     plot(t_array,delta)
%     xlabel('Time')
%     ylabel('delta')


end

%Direct differentiation

    function sensitivity = direct_diff(k1,k2,t_total,h,scale_factor )
        y0 = 0.0;
        alphaMax = 0;
        t = 0;
        obj_function = 0;

        m = 1;

        
        
        y = zeros(2,1);
        Dy = zeros(2,1);
        Dq = 0;
        Dp = 1;
        DalphaMax = 0;

        N_steps = length(0:h:t_total);
        delta = zeros(1,N_steps);
        t_array = zeros(1,N_steps);
        force_array = zeros(1,N_steps);
        k = 1;
        for i=0:h:t_total-h

            % Midpoint method
        
            % First evaluation            
            delta_y2 = y(2);
            dFdY = TangentdY(delta_y2,alphaMax,k1,k2,m);
            dFdY_global = [0 , dFdY; 1, 0];
            dFdAlphaMax = TangentdAlphaMax(delta_y2,alphaMax,k1,k2,m);
            dFdAlphaMax_global = [dFdAlphaMax; 0];
            dFdP = TangentdP(delta_y2,alphaMax,k1,k2,m);
            dFdP_global = [dFdP ; 0];
            dRdY = ObjdY(delta_y2,alphaMax,k1,k2,m,scale_factor );
            dRdY_global = [0; dRdY];
            
            [force, plast,int_force] = RHS(y(2),alphaMax,k1,k2,m,t);
            
            %Update the variables
            %obj_function = obj_function + 1/2*h*y(2)/scale_factor;
            y(2) = y(2) + 1/2*h*y(1);
            y(1) = y(1) + 1/2*h*force;
        
            
            
            %Dq = Dq + 1/2*h* dRdY_global'*Dy;
            Dy = Dy + 1/2*h*( dFdY_global*Dy + dFdAlphaMax_global*DalphaMax + dFdP_global * Dp);
%             disp('Tiempo')
%             disp(i)
%             disp('delta')
%             disp(delta_y2)
%             disp('AlphaMax')
%             disp(alphaMax)
            
            
            % Value we'll feed to the next stage
            t = t + 1/2*h; 
        
            % Call again the RHS
            [force, plast,int_force] = RHS(y(2),alphaMax,k1,k2,m,t);
            
            % Second Evaluation
            
            delta_y2 = y(2);
            dFdY = TangentdY(delta_y2,alphaMax,k1,k2,m);
            dFdY_global = [0 , dFdY; 1, 0];
            dFdAlphaMax = TangentdAlphaMax(delta_y2,alphaMax,k1,k2,m);
            dFdAlphaMax_global = [dFdAlphaMax; 0];
            dFdP = TangentdP(delta_y2,alphaMax,k1,k2,m);
            dFdP_global = [dFdP ; 0];
            dRdY = ObjdY(delta_y2,alphaMax,k1,k2,m,scale_factor );
            dRdY_global = [0; dRdY];
        
            


            % Update the variable
            obj_function = obj_function + h*y(2)/scale_factor;
            y(2) = y(2) + h*y(1);
            y(1) = y(1) + h*force;

            Dq = Dq + h* dRdY_global'*Dy;
            Dy = Dy + h*( dFdY_global*Dy + dFdAlphaMax_global*DalphaMax + dFdP_global * Dp);
            
            % Update Plastic variable
            % Call again the RHS
            %[force, plast,int_force] = RHS(y(2),alphaMax,k1,k2,m,t);
            if plast == 1
                dHdY = 1;
                dHdalpha = 0;
                alphaMax = y(2);
            else
                dHdY = 0;
                dHdalpha = 1;
            end
            
            
            DalphaMax = dHdalpha*DalphaMax + dHdY*Dy(2);


            t = t + 1/2*h; 
            delta(k) = y(2);
            t_array(k) = t;
            force_array(k) = int_force;
            
            k = k + 1;

            

        end
        
        obj_function
        
        
        sensitivity_ultimo = Dy(2);
        
%         figure(3)
%         plot(delta,force_array)
%     
%         figure(4)
%         plot(t_array,force_array)
        disp('Sensitivity')
        disp(Dq)
        
%         disp('Sensitivity ultimo')
%         disp(sensitivity_ultimo)
    end

    function sensitivity = ...
            adjoint_method (k1,k2,t_total,h, delta_history, alphaMaxhistory,scale_factor,delta_intermedio)
        
        y = zeros(2,1);
        lambda = zeros(2,1);
        beta = 0;
        gamma = 0;
        m = 1;

        N_steps = length(0:h:t_total);
        delta = zeros(1,N_steps);
        t_array = zeros(1,N_steps);
        gamma_array = zeros(1,N_steps);
        k = 0;
        for i=t_total-h:-h:0

            delta_y2 = delta_intermedio(N_steps - k - 1);
            alphaMax = alphaMaxhistory(N_steps - k - 1);
            dFdY = TangentdY(delta_y2,alphaMax,k1,k2,m);
            dFdY_global = [0 , dFdY; 1, 0];
            dFdAlphaMax = TangentdAlphaMax(delta_y2,alphaMax,k1,k2,m);
            dFdAlphaMax_global = [dFdAlphaMax; 0];
            dFdP = TangentdP(delta_y2,alphaMax,k1,k2,m);
            dFdP_global = [dFdP ; 0];
            dRdY = ObjdY(delta_y2,alphaMax,k1,k2,m,scale_factor);
            dRdY_global = [0; dRdY];
            
            if dFdAlphaMax ~= 0
                a = 1;
            end
            % Calculate U_2 
            U_2 = h * dFdY_global * lambda + h*dRdY_global;
            
            V_2 = h * dFdP_global' * lambda;
            
            W_2 = h * dFdAlphaMax_global' * lambda;
            
%             disp('Tiempo')
%             disp(t_total-i)
%             disp('delta')
%             disp(delta_y2)
%             disp('AlphaMax')
%             disp(alphaMax)


            delta_y2 = delta_history(N_steps - k - 1);
            alphaMax = alphaMaxhistory(N_steps - k - 1);
            dFdY = TangentdY(delta_y2,alphaMax,k1,k2,m);
            dFdY_global = [0 , dFdY; 1, 0];
            dFdAlphaMax = TangentdAlphaMax(delta_y2,alphaMax,k1,k2,m);
            dFdAlphaMax_global = [dFdAlphaMax; 0];
            dFdP = TangentdP(delta_y2,alphaMax,k1,k2,m);
            dFdP_global = [dFdP ; 0];
            dRdY = ObjdY(delta_y2,alphaMax,k1,k2,m,scale_factor);
            dRdY_global = [0; dRdY];

            if dFdAlphaMax ~= 0
                a = 1;
            end
            U_1 = h * dFdY_global * U_2 * 1/2;
            
            V_1 = h * dFdP_global'* U_2 * 1/2;
            
            W_1 = h * dFdAlphaMax_global'* U_2 *  1/2;
            
            lambda = lambda + U_2 + U_1;
            
            beta = beta + V_2 + V_1;
            
            % Call state equation to check for update
            delta = delta_history(N_steps - k);
            dHdalpha = state_eq_alphaMax(delta,alphaMax,k1,k2,m);
            
            gamma = gamma * dHdalpha + W_1 + W_2;
            if (i > h)
                delta_anterior = delta_history(N_steps - k - 1);
                alphaMax_anterior = alphaMaxhistory(N_steps - k - 2);
                
%                 disp('delta anterior')
%                 disp(delta_anterior)
%                 disp('AlphaMax anterior')
%                 disp(alphaMax_anterior)
                dHdY = state_eq_Y(delta_anterior,alphaMax_anterior,k1,k2,m);
            else
                dHdY = 0;
            end
            
            dHdY_global = [0; dHdY];

            lambda = lambda + dHdY_global*gamma;
            


            k = k + 1;
        end    
        
        
        disp('Sensitivity adjoint')
        disp(beta)
        sensitivity = beta;
        
%         disp('Sensitivity ultimo adjoint')
%         disp(sensitivity)
        
        
        
    end

    function dFdY = TangentdY(delta,alphaMax,k1,k2,m)
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1*delta;
        

        if ElasticForce >= PlasticForce
           dFdY = k1;
           alphaMax = delta;
        else
           dFdY = k2;             
        end      
        
        dFdY = - dFdY* 1/m;
        
    end

    function dFdAlphaMax = TangentdAlphaMax(delta,alphaMax,k1,k2,m)
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1*delta;
        
        if ElasticForce >= PlasticForce
           dFdAlphaMax = 0;
           alphaMax = delta;
        else
           dFdAlphaMax = k2*(- (k2 - k1)/ k2);             
        end      
        
        dFdAlphaMax = - dFdAlphaMax* 1/m;
        
    end

    function dFdP = TangentdP(delta,alphaMax,k1,k2,m)
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1*delta;
        
        if ElasticForce >= PlasticForce
           dFdP = delta;
           alphaMax = delta;
        else
           dFdP = alphaMax;             
        end      
        
        dFdP = - dFdP* 1/m;
        
    end

    function dRdY = ObjdY(delta,alphaMax,k1,k2,m,scale_factor)
        
        dRdY = 1/scale_factor;
        
    end

    function dHdalpha = state_eq_alphaMax(delta,alphaMax,k1,k2,m)
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1*delta;
        
        if ElasticForce >= PlasticForce
           dHdalpha = 0;
           alphaMax = delta;
        else
           dHdalpha = 1;             
        end      
        
        
    end

    function dHdY = state_eq_Y(delta,alphaMax,k1,k2,m)
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1*delta;
        
        if ElasticForce >= PlasticForce
           dHdY = 1;
           alphaMax = delta;
        else
           dHdY = 0;             
        end      
        
        
    end


    function [force,plast,int_force] = RHS(delta,alphaMax,k1,k2,m,t)

    %Call external force
    ext_force = ForceExt(t);

    %Reinitialize to zero, write 1 if plastic regime
    plast = 0;
    %Call internal force
    [int_force, plast] = ForceInt(delta,alphaMax,k1,k2);

    force = 1/m * (ext_force - int_force);

    end

    function ext_force = ForceExt(t)
        
        
        if t < 20
            ext_force = 2*t;
        else
           ext_force = 80 -t*2;
        end
        
    end

    function [int_force,plast] = ForceInt(delta,alphaMax,k1,k2)
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1*delta;
        
        tol = 1e-15;
        if abs(ElasticForce - PlasticForce) < tol
            int_force = PlasticForce;
            plast = 1;            
        elseif ElasticForce >= PlasticForce
            int_force = PlasticForce;
            plast = 1;
        else
           int_force = ElasticForce; 
           plast = 0;
            
        end
    end
end