function dummy

global tol;

tol = 1e-10;
n_elem = 2;
initial_design = 0.5;

design_beads = [1:n_elem+1];
design_array = zeros(1,n_elem+1);

fd = 1;

for l=1:n_elem+1
    
    if any(l == design_beads)
        design_array(l) = 1;
    end
    
end

design_array = design_array == 1;

k1 = zeros(n_elem,1);

densities = ones(1,n_elem+1);
densities(design_array) = initial_design;
  
densities_matrix = [densities(1:end-1); densities(2:end)];
k1 = mean(densities_matrix,1);

k2 = 4;
h = 0.2;
scale_factor = 1;

format long


t_total = 6.206;
%t_total = 4.206;
t_total = 250
[obj_function_1,obj_function_ultimo,y_history,alphaMaxhistory ] ...
    = primal_problem(k1,k2,t_total,h,scale_factor,n_elem, densities,design_array);

obj_function_1
obj_function_ultimo

[sensitivity_diff,sensitivity_ultimo] = direct_diff(k1,k2,t_total,h,scale_factor,n_elem,densities, design_array,design_beads);
% % 
% % 
[sensitivity_adjoint,sensitivity_ultimo_adjoint] = adjoint_method (k1,k2,t_total,h, y_history, alphaMaxhistory,scale_factor,n_elem,densities,design_array,design_beads);
% % 
% disp('Function analitica')
% phi = 2 / k1^2 * cos(sqrt(k1) * t_total) + t_total^2 /k1 - 2 / k1^2;
% disp(ph   i)
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

if fd == 1
    sensitivity_FD = zeros(1,sum(design_array));
    sensitivity_FD_ultimo = zeros(1,sum(design_array));
    contador = 1;
    for q = design_beads
        dt = 1e-7;
        densities(q)  = densities(q)  + dt;

        densities_matrix = [densities(1:end-1); densities(2:end)];
        k1 = mean(densities_matrix,1);

        [obj_function_2,obj_function_ultimo_2,y_history, alphaMaxhistory] ...
            = primal_problem(k1,k2,t_total,h,scale_factor,n_elem,densities,design_array);

        densities(q)  = densities(q)  - 2*dt;

        densities_matrix = [densities(1:end-1); densities(2:end)];
        k1 = mean(densities_matrix,1);

        [obj_function_3, obj_function_ultimo_3, y_history,alphaMaxhistory]...
            = primal_problem(k1,k2,t_total,h,scale_factor,n_elem,densities,design_array);


        densities(q)  = densities(q)  + dt;
        densities_matrix = [densities(1:end-1); densities(2:end)];
        k1 = mean(densities_matrix,1);
        sensitivity_FD(contador) = (obj_function_2 - obj_function_3)/(2*dt);
        sensitivity_FD_ultimo(contador) = (obj_function_ultimo_2 - obj_function_ultimo_3)/(2*dt);
        contador = contador + 1;
    end

    disp('Sensitivity Comparison')
    disp('FD               Adjoint                 Direct Diff')
    disp([sensitivity_FD', sensitivity_adjoint, sensitivity_diff] )

%     disp('Sensitivity ultimo Comparison')
%     disp('FD               Adjoint                 Direct Diff')
%     disp([sensitivity_FD_ultimo', sensitivity_ultimo_adjoint, sensitivity_ultimo] )
end

%sensitivity_FD_ultimo = (obj_function_ultimo_2 - obj_function_ultimo_3)/(2*dt)




function [obj_function, obj_function_ultimo, y_history, alphaMaxhistory] = ...
        primal_problem(k1,k2,t_total,h,scale_factor,n_elem,densities,design_array)
    
    alphaMax = zeros(n_elem,1);
    t = 0;
    obj_function = 0;

    
    R = 4.7629999999999999;
    
    mass_array =  densities*4.0/3.0*pi*R^3*8.5;
    m = 4.0/3.0*pi*R^3*8.5;


    y = zeros(2*n_elem,1);

    N_steps = length(0:h:t_total);
    delta = zeros(n_elem,N_steps);
    alphaMaxhistory = zeros(n_elem,N_steps);
    
    
    t_array = zeros(1,N_steps);
    force_array = zeros(n_elem,N_steps);
    k = 1;
    plast = 0;
    
    plast_elem = zeros(n_elem,1);
    int_force_elem = zeros(n_elem,1);
    y_elem_array = zeros(n_elem,1);
    int_force = 0;
    
    
    dof_elements = [1:2:(n_elem*2-1); 2:2:(n_elem*2); 3:2:(n_elem*2 + 1); 4:2:(n_elem*2 + 2)];

    y_history = zeros(2*n_elem,N_steps);
    for i=0:h:t_total
        

        y_history(:,k) = y;
          
        y_asemmbly = zeros(2*n_elem,1);
        for j = 1:n_elem

            delta(j,k) = y_elem_array(j);
            if plast_elem(j) == 1
                alphaMax(j) = y_elem_array(j);
            end
            alphaMaxhistory(j,k) = alphaMax(j);
            force_array(j,k) = int_force_elem(j);
            
            if j ~= n_elem
                y_elem = -1.0*diff(y(dof_elements(2:2:end,j)));
            else
                y_elem = y(dof_elements(2,j));
            end
            
%             if i == 3
%                 disp('Primal at 3')
%                 disp(j)
%                 disp(y_elem)
%                 disp(alphaMax(j))
%             end
            

            [force, plast,int_force] = RHS(y_elem,alphaMax(j),k1,k2,mass_array,t,j);

            if j ~= n_elem
                y_asemmbly(dof_elements(2:2:end,j)) = h*y(dof_elements(1:2:end,j));
                y_asemmbly(dof_elements(1:2:end,j)) = y_asemmbly(dof_elements(1:2:end,j)) + h*force;
            else
                y_asemmbly(dof_elements(2,j)) = h*y(dof_elements(1,j));
                y_asemmbly(dof_elements(1,j)) = y_asemmbly(dof_elements(1,j)) + h*force(1);                  
            end          
            y_elem_array(j) = y_elem;


            
%             if j == 2
%                 disp(y_elem)
%                 disp(int_force)
%             end


            if j == n_elem
                obj_function = obj_function + h*y(dof_elements(2,j))/scale_factor;
            end
        end
        y = y + y_asemmbly;
        
        for j=1:n_elem
            if j ~= n_elem
                y_elem = -1.0*diff(y(dof_elements(2:2:end,j)));
            else
                y_elem = y(dof_elements(2,j));
            end
            
            y_elem_array(j) = y_elem;
            
            %Update alphaMax
            [force, plast,int_force] = RHS(y_elem,alphaMax(j),k1,k2,mass_array,t,j);
            plast_elem(j) = plast;
            int_force_elem(j) = int_force;
            
        end
        
        t = t + h; 
        t_array(k) = t;
        k = k + 1;
    end

    % Last update
    delta(:,k) = y_elem_array ;
    
    for j = 1:n_elem
        if plast_elem(j) == 1
            alphaMax(j) = y_elem_array(j);
        end
        force_array(j,k) = int_force_elem(j);
    end
    alphaMaxhistory(:,k) = alphaMax;
    
    y_history(:,k) = y;
    
%     ll = 1;
%     
%     figure('OuterPosition',[0 0  1600 1200])
%     writerObj = VideoWriter('peaks.avi');
%     open(writerObj);
%     init_coord = [1:1:n_elem+1];
%     for i=0:h:t_total-h
%        cla
%        position = init_coord + 1e+1*y_history(2:2:end,ll);
%        scatter(position, ones(n_elem+1,1));
%        hold on
%        
%        if mod(ll,2)==0
%             frame = getframe;
%             writeVideo(writerObj,frame);
%         end
%        ll = ll + 1;
%     end
%     
%     close(writerObj);
    
    for i=1:n_elem
        plot(delta(i,1:end),force_array(i,1:end))
        xlabel('delta')
        ylabel('force')
        hold all
    end
    

%     
%     figure(2)
%     plot(t_array,alphaMaxhistory)
    
    
    obj_function_ultimo = y(2*(n_elem-1));
    
    
%     figure
%     plot(t_array,delta)
%     xlabel('Time')
%     ylabel('delta')


end

%Direct differentiation

    function [sensitivity,sensitivity_ultimo] = direct_diff(k1,k2,t_total,h,scale_factor,n_elem , densities, design_array,design_beads)

        

        R = 4.7629999999999999;
        rho = 8.5;
        mass_array =  densities*4.0/3.0*pi*R^3*8.5;
        deriv_mass = 4.0/3.0*pi*R^3*rho;
        
        n_variables = sum(design_array);
        sensitivity = zeros(n_variables,1);
        sensitivity_ultimo = zeros(n_variables,1);
        
        contador_diff = 1;
        for variable = design_beads
        
            t = 0;
            alphaMax = zeros(n_elem,1);
            obj_function = 0;
            y = zeros(2*n_elem,1);
            Dy = zeros(2*n_elem,1);
            Dq = 0;
            Dp = 1;
            DalphaMax =  zeros(n_elem,1);

            y_elem_array = zeros(n_elem,1);


            N_steps = length(0:h:t_total);
            delta = zeros(n_elem,N_steps);
            t_array = zeros(1,N_steps);
            force_array = zeros(n_elem,N_steps);
            int_force_elem = zeros(n_elem,1);
            k = 1;
            plast = 0;

            dof_elements = [1:2:(n_elem*2-1); 2:2:(n_elem*2); 3:2:(n_elem*2 + 1); 4:2:(n_elem*2 + 2)];

            plast_elem = zeros(n_elem,1);
            for i=0:h:t_total


                y_asemmbly = zeros(2*n_elem,1);
                Dy_asemmbly = zeros(2*n_elem,1);
                for j=1:n_elem
                    delta(j,k) = y_elem_array(j);
                    if plast_elem(j) == 1
                        alphaMax(j) = y_elem_array(j);
                    end
                    force_array(j,k) = int_force_elem(j);


                    if j ~= n_elem
                        y_elem = -1.0*diff(y(dof_elements(2:2:end,j)));
                    else
                        y_elem = y(dof_elements(2,j));
                    end

                    dFdY = TangentdY(y_elem,alphaMax(j),k1,k2,mass_array,j);
                    dFdY_global = [0 , dFdY(1,1), 0 dFdY(1,2);...
                                   1, 0, 0, 0;...
                                   0, dFdY(2,1), 0 , dFdY(2,2); ...
                                   0, 0, 1, 0];
                    dFdAlphaMax = TangentdAlphaMax(y_elem,alphaMax(j),k1,k2,mass_array,j);
                    dFdAlphaMax_global = [dFdAlphaMax(1); 0; dFdAlphaMax(2); 0];
                    dFdP = TangentdP(y_elem,alphaMax(j),k1,k2,mass_array,j,n_elem,deriv_mass,design_array,variable,t);
                    dFdP_global = [dFdP(1); 0; dFdP(2); 0];


        %             disp('Tiempo')
        %             disp(i)
        %             disp('delta')
        %             disp(delta_y2)
        %             disp('AlphaMax')
        %             disp(alphaMax)

                    [force, plast,int_force] = RHS(y_elem,alphaMax(j),k1,k2,mass_array,t,j);



                    if j ~= n_elem
                        y_asemmbly(dof_elements(2:2:end,j)) =  h*y(dof_elements(1:2:end,j));
                        y_asemmbly(dof_elements(1:2:end,j)) = y_asemmbly(dof_elements(1:2:end,j)) ...
                                                            + h*force;

                        Dy_asemmbly(dof_elements(1:2:end,j)) = Dy_asemmbly(dof_elements(1:2:end,j)) ...
                                        + h*( dFdY_global(1:2:end,2:2:end)*Dy(dof_elements(2:2:end,j)) ...
                                        + dFdAlphaMax_global(1:2:end)*DalphaMax(j) ...
                                        + dFdP_global(1:2:end) * Dp);                                

                        Dy_asemmbly(dof_elements(2:2:end,j)) = h*dFdY_global(2:2:end,1:2:end)*Dy(dof_elements(1:2:end,j));       

                        Dy_asemmbly(dof_elements(2:2:end,j)) = Dy_asemmbly(dof_elements(2:2:end,j)) ...
                                        + h*(dFdAlphaMax_global(2:2:end)*DalphaMax(j) ...
                                                + dFdP_global(2:2:end) * Dp);
                    else
                        y_asemmbly(dof_elements(2,j)) = h*y(dof_elements(1,j));
                        y_asemmbly(dof_elements(1,j)) = y_asemmbly(dof_elements(1,j))...
                                                        + h*force(1);     

                        aux_vector = dFdY_global(1:2,1:2)*Dy(dof_elements(1:2,j));

                        Dy_asemmbly(dof_elements(1,j)) = Dy_asemmbly(dof_elements(1,j)) ...
                                        + h*( aux_vector(1) + dFdAlphaMax_global(1)*DalphaMax(j) ...
                                        + dFdP_global(1) * Dp);

                        Dy_asemmbly(dof_elements(2,j)) = h*aux_vector(2);            
                        Dy_asemmbly(dof_elements(2,j)) = Dy_asemmbly(dof_elements(2,j)) ...
                                        + h*(dFdAlphaMax_global(2)*DalphaMax(j) ...
                                        + dFdP_global(2) * Dp);
                    end          
                    y_elem_array(j) = y_elem;



                    [force, plast,int_force] = RHS(y_elem,alphaMax(j),k1,k2,mass_array,t,j);


                    int_force_elem(j) = int_force;


                    plast_elem(j) = plast;
                    if j == n_elem
                        obj_function = obj_function + h*y(dof_elements(2,j))/scale_factor;
                        dRdY = ObjdY(y_elem,alphaMax(j),k1,k2,scale_factor );
                        dRdY_global = [0; dRdY];
                        Dq = Dq + h* dRdY_global'*Dy(dof_elements(1:2,j));
                    end

                end

                y = y +  y_asemmbly;

                Dy = Dy + Dy_asemmbly;

                for j=1:n_elem
                    if j ~= n_elem
                        y_elem = -1.0*diff(y(dof_elements(2:2:end,j)));
                    else
                        y_elem = y(dof_elements(2,j));
                    end

                    y_elem_array(j) = y_elem;

                    %Update alphaMax
                    [force, plast,int_force] = RHS(y_elem,alphaMax(j),k1,k2,mass_array,t,j);
                    plast_elem(j) = plast;
                    int_force_elem(j) = int_force;

                end

                for j = 1:n_elem


                    if j ~= n_elem
                        if plast_elem(j) == 1
                            dHdY = [1, -1];
                            dHdalpha = 0;
                        else
                            dHdY = [0, 0];
                            dHdalpha = 1;                
                        end

                        DalphaMax(j) = dHdalpha*DalphaMax(j) + dHdY*Dy(dof_elements(2:2:end,j)); 
                    else                    
                        if plast_elem(j) == 1
                            dHdY = 1;
                            dHdalpha = 0;
                        else
                            dHdY = 0;
                            dHdalpha = 1;                
                        end

                        DalphaMax(j) = dHdalpha*DalphaMax(j) + dHdY*Dy(dof_elements(2,j));                          
                    end

                end

                k = k + 1;

                t = t + h;

            end

            % Last update
            delta(:,k) = y_elem_array ;

            for j = 1:n_elem
                if plast_elem(j) == 1
                    alphaMax(j) = y_elem_array(j);
                end
                force_array(j,k) = int_force_elem(j);
            end

    %         plot(delta(1,1:end),force_array(1,1:end))
    %         xlabel('delta')
    %         ylabel('force')
    %         hold on
    %          plot(delta(2,1:end),force_array(2,1:end),'r')
            sensitivity_ultimo(contador_diff) = Dy(2*(n_elem-1));

    %         figure(3)
    %         plot(delta,force_array)
    %     
    %         figure(4)
    %         plot(t_array,force_array)
            obj_function
            sensitivity(contador_diff) = Dq;
            contador_diff = contador_diff + 1;
        end
    end

    function [sensitivity,sensitivity_ultimo] = adjoint_method (k1,k2,t_total,h, y_history, ...
                        alphaMaxhistory,scale_factor,n_elem,densities,design_array,design_beads )
        

        n_variable = sum(design_array);
        
        sensitivity_ultimo = zeros(n_variable,1);
        alphaMax = zeros(n_elem,1);
        
        R = 4.7629999999999999;
        rho = 8.5;
        mass_array =  densities*4.0/3.0*pi*R^3*8.5;
        deriv_mass = 4.0/3.0*pi*R^3*rho;

        y = zeros(2*n_elem,1);

        N_steps = length(0:h:t_total);
        
        lambda = zeros(2*n_elem,1);
        
        lambda_array = zeros(2*n_elem,N_steps);
        beta = zeros(n_variable,1);
        
        beta_array = zeros(n_variable,N_steps);
        gamma = zeros(n_elem,1);

        plast_elem = zeros(n_elem,1);
        int_force_elem = zeros(n_elem,1);
        y_elem_array = zeros(n_elem,1);
        int_force = 0;

        dof_elements = [1:2:(n_elem*2-1); 2:2:(n_elem*2); 3:2:(n_elem*2 + 1); 4:2:(n_elem*2 + 2)];

        gamma_array = zeros(n_elem,N_steps);
        k = N_steps ;
        
        t = t_total + h;

        for k=N_steps:-1:1
            
            t = t - h;

            lambda_assembly = zeros(2*n_elem,1);
            for j=1:n_elem
                
                % Grab history values
                alphaMax = alphaMaxhistory(j,k);                
                if j ~= n_elem
                    y_elem = -1.0*diff(y_history(dof_elements(2:2:end,j),k));
                else
                    y_elem = y_history(dof_elements(2,j),k);
                end
                
%                 if i == 3
%                     disp('Adjoint at 3')
%                     disp(j)
%                     disp(y_elem)
%                     disp(alphaMax)
%                 end
                dFdY = TangentdY(y_elem,alphaMax,k1,k2,mass_array,j);
                dFdY_global = [0 , dFdY(1,1), 0 dFdY(1,2);...
                               1, 0, 0, 0;...
                               0, dFdY(2,1), 0 , dFdY(2,2); ...
                               0, 0, 1, 0];


                
                if j ~= n_elem
                    trans_dFdY_global = dFdY_global';
                    lambda_assembly(dof_elements(1:2:end,j)) = ...
                                h*(trans_dFdY_global(1:2:end,2:2:end) * lambda(dof_elements(2:2:end,j)));
                    
                    lambda_assembly(dof_elements(2:2:end,j)) = ...
                                lambda_assembly(dof_elements(2:2:end,j))...
                                + h*(trans_dFdY_global(2:2:end,1:2:end) * lambda(dof_elements(1:2:end,j)));

                else
                    dRdY = ObjdY(y_elem,alphaMax,k1,k2,scale_factor );
                    dRdY_global = [0; dRdY];
                    
                    aux_mat = dFdY_global';

                    lambda_assembly(dof_elements(1,j)) = ...
                                h*(aux_mat(1,2) * lambda(dof_elements(2,j)) + dRdY_global(1)); 
                    lambda_assembly(dof_elements(2,j)) = ...
                                lambda_assembly(dof_elements(2,j))...
                                + h*(aux_mat(2,1) * lambda(dof_elements(1,j)) + dRdY_global(2));               
                end                  

            end   
            
            
            contador_adjoint = 1;
            for variable = design_beads
                for j=1:n_elem

                    % Grab history values
                    alphaMax = alphaMaxhistory(j,k);                
                    if j ~= n_elem
                        y_elem = -1.0*diff(y_history(dof_elements(2:2:end,j),k));
                    else
                        y_elem = y_history(dof_elements(2,j),k);
                    end
        
                    dFdP = TangentdP(y_elem,alphaMax,k1,k2,mass_array,j,n_elem,deriv_mass,design_array,variable,t);
                    dFdP_global = [dFdP(1); 0; dFdP(2); 0];

                    if j ~= n_elem
                        beta(contador_adjoint) = beta(contador_adjoint) + h*dFdP_global'*lambda(dof_elements(:,j));
                    else
                        dFdP_global_aux = dFdP_global(1:2,:);
                        beta(contador_adjoint) = beta(contador_adjoint) + h*dFdP_global_aux'*lambda(dof_elements(1:2,j));
                    end


                end
                contador_adjoint = contador_adjoint + 1;
            end
            
            % State equation evaluated at the displacement calculated at
            % the current time step  
            for j=1:n_elem
 
                alphaMax = alphaMaxhistory(j,k);
                if j ~= n_elem
                    y_elem = -1.0*diff(y_history(dof_elements(2:2:end,j),k + 1));
                else
                    y_elem = y_history(dof_elements(2,j),k + 1);
                end
                dHdalpha = state_eq_alphaMax(y_elem,alphaMax,k1,k2);    

                
                gamma(j) = gamma(j) * dHdalpha;
                
            end
            
            
            for j=1:n_elem       
                
                % Grab history values
                alphaMax = alphaMaxhistory(j,k);                
                if j ~= n_elem
                    y_elem = -1.0*diff(y_history(dof_elements(2:2:end,j),k));
                else
                    y_elem = y_history(dof_elements(2,j),k);
                end
                if j == 1 && N_steps - k + 1 == 611 
                    a= 1;
                end
                
                dFdAlphaMax = TangentdAlphaMax(y_elem,alphaMax,k1,k2,mass_array,j);
                
                dFdAlphaMax_global = [dFdAlphaMax(1); 0; dFdAlphaMax(2); 0];
                

                if j ~= n_elem
                    gamma(j) = gamma(j)...
                                    + h * ( dFdAlphaMax_global' * lambda(dof_elements(:,j)));
                else
                    dFdAlphaMax_global = dFdAlphaMax_global(1:2);
                    gamma(j) = gamma(j) ...
                                    + h * ( dFdAlphaMax_global' * lambda(dof_elements(1:2,j)));
                end                
            end
            
            
            lambda(:) = lambda(:) + lambda_assembly;
            
            
            % Call state equation to calculate the gamma
            % contribution to lambda
            if (k > 1)             
                
                for j=1:n_elem 
                    alphaMax = alphaMaxhistory(j,k - 1);
                    if j ~= n_elem
                        y_elem = -1.0*diff(y_history(dof_elements(2:2:end,j),k));
                    else
                        y_elem = y_history(dof_elements(2,j),k);
                    end            
                    dHdY = state_eq_Y(y_elem,alphaMax,k1,k2);



                    if j ~= n_elem
                        lambda(dof_elements(2:2:end,j)) = ...
                            lambda(dof_elements(2:2:end,j)) + dHdY*[1;-1]*gamma(j);
                    else
                        lambda(dof_elements(2,j)) = ...
                            lambda(dof_elements(2,j)) + dHdY*gamma(j);
                    end           
                end
            else
                dHdY = 0;
            end
            
            gamma_array(:,N_steps - k + 1) = gamma;
            
            lambda_array(:,N_steps - k + 1) = lambda;
            
            beta_array(:,N_steps - k + 1) = beta;

        end
        
    sensitivity = beta;
        
        
        
        
        
    end

    function dFdY = TangentdY(delta,alphaMax,k1_array,k2,mass_array,j)
        
        k1_elem = k1_array(j);
        alphaP = alphaMax*(k2 - k1_elem)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1_elem*delta;

        tol = 1e-10;
        if delta >= alphaP
            
            if ElasticForce < PlasticForce - tol
                dFdY = k2;       
            else
               dFdY = k1_elem;
               alphaMax = delta;         
            end      
        else
            dFdY = 0;
        end
        
        dFdY = [dFdY, -dFdY;-dFdY, dFdY];
        
        dFdY(1,:) = - dFdY(1,:)* 1/mass_array(j);
        dFdY(2,:) = - dFdY(2,:)* 1/mass_array(j+1);
            
        
    end

    function dFdAlphaMax = TangentdAlphaMax(delta,alphaMax,k1_array,k2,mass_array,j)
        
        k1_elem = k1_array(j);
        alphaP = alphaMax*(k2 - k1_elem)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1_elem*delta;
        
        if delta >= alphaP
            
            if ElasticForce < PlasticForce - tol
                dFdAlphaMax = k2*(- (k2 - k1_elem)/ k2);             
            else
                dFdAlphaMax = 0;
                alphaMax = delta;       
            end      
        else
            dFdAlphaMax = 0;
        end
        
        dFdAlphaMax = - dFdAlphaMax;
        
        dFdAlphaMax = [dFdAlphaMax* 1/mass_array(j); -dFdAlphaMax* 1/mass_array(j+1)];
        
    end

    function dFdP_array = TangentdP(delta,alphaMax,k1_array,k2,mass_array,j,n_elem,deriv_mass,design_array, variable,t)
        
        k1_elem = k1_array(j);
        
        alphaP = alphaMax*(k2 - k1_elem)/ k2;
        
        DalphaP = alphaMax/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1_elem*delta;
        
        dFdP_array = zeros(2,1);
        
        dFdP = 0;
        
        ext_force = ForceExt(t);
        if j ~= 1
            ext_force = 0;
        end
        
        if design_array(j) && j == variable
            if delta >= alphaP
                if ElasticForce < PlasticForce - tol
                   F = ElasticForce;
                   dFdP = alphaMax*1/2;  
                else
                   dFdP = delta*1/2;
                   F = PlasticForce;
                   alphaMax = delta;     
                end      
            else
                dFdP = 0;
                F = 0;
            end
            
            dFdP_array(1) = (-1.0*(ext_force-F)*deriv_mass/mass_array(j)^2-dFdP* 1/mass_array(j));
            
            dFdP_array(2) = dFdP* 1/mass_array(j+1);
        
        end
        
        if design_array(j+1) && j + 1 == variable
        
            if delta >= alphaP
                if ElasticForce < PlasticForce - tol
                   F = ElasticForce;
                   dFdP = alphaMax*1/2;  
                else
                   dFdP = delta*1/2;
                   F = PlasticForce;
                   alphaMax = delta;       
                end      
            else
                dFdP = 0;
                F = 0;
            end
            
            dFdP_array(1) = -dFdP* 1/mass_array(j);
            
            dFdP_array(2) = (-1.0*(F)*deriv_mass/mass_array(j+1)^2+dFdP* 1/mass_array(j+1));
            
            
        
        end
        
         
    end

    function dRdY = ObjdY(delta,alphaMax,k1,k2,scale_factor)
        
        dRdY = 1/scale_factor;
        
    end

    function dHdalpha = state_eq_alphaMax(delta,alphaMax,k1,k2)
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1*delta;
           
        if ElasticForce < PlasticForce
           dHdalpha = 1;      
        else
           dHdalpha = 0;
           alphaMax = delta;           
        end      
        
        
    end

    function dHdY = state_eq_Y(delta,alphaMax,k1,k2)
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        ElasticForce = k2*(delta - alphaP);
        
        PlasticForce = k1*delta;
        
        
        if ElasticForce < PlasticForce
           dHdY = 0;
        else
           dHdY = 1;
           alphaMax = delta;           
        end      
        
        
    end


    function [force,plast,int_force] = RHS(delta,alphaMax,k1_array,k2,mass_array,t,j)

        %Call external force
        ext_force = ForceExt(t);
        if j ~= 1
            ext_force = 0;
        end

        %Reinitialize to zero, write 1 if plastic regime
        plast = 0;
        %Call internal force
        k1_elem = k1_array(j);
        [int_force, plast] = ForceInt(delta,alphaMax,k1_elem,k2);


            
       
        force = [1/mass_array(j) * (ext_force - int_force); 1/mass_array(j+1) * (int_force)];

    end

    function ext_force = ForceExt(t)
        
        
        if t < 40
           ext_force = 10*t;
        elseif t >= 40 && t < 80
           ext_force = 800 -t*10;
        else 
           ext_force = 0;
        end
                
        
    end

    function [int_force,plast] = ForceInt(delta,alphaMax,k1,k2)
        
        
        alphaP = alphaMax*(k2 - k1)/ k2;
        
        if delta >= alphaP
        
            ElasticForce = k2*delta - alphaMax*(k2 - k1);

            PlasticForce = k1*delta;

            if ElasticForce < PlasticForce - tol
               int_force = ElasticForce; 
               plast = 0;
            else
               int_force = PlasticForce;
               plast = 1;
            end
        else
            int_force = 0;
            plast = 0;
        end
    end
end