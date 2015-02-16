function [delta, alphaMax, force,...
    new_delta,new_force, new_alphaMax, ...
    alphaMaxUnloading,delta_unloaded,...
    delta_plastic_regime, force_plastic_regime,...
    delta_elastic_regime, force_elastic_regime, alphaMax_elastic_regime,...
    smoothed_delta, smoothed_force] = validate_data(plotting,data,j)

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
R_from = 0.0020*1e3;
%Combined radius
RStar = (R_to*R_from)/(R_to+R_from);

%Yield stress
sigmayielding = 0.55; %Gpa

% Yield parameters
alphaY = (pi/2)^2* (sigmayielding*1.6/EStar)^2 * RStar;


M = data;    

if plotting == 1
    plot(M(:,1),M(:,2))
    hold on
end
smooth_parameter = 50;
yy2 = smooth(M(:,2),smooth_parameter);
xx2 = smooth(M(:,1),smooth_parameter); 


if plotting == 1
    plot(xx2,yy2, 'r*')

    xlabel('Strain')
    ylabel('Force')
    legend('Experimental data','Smoothed experimental data')
end
smoothed_delta = xx2;
smoothed_force = yy2;


threshold_parameter = 0;
% Find local minima for strain
[strainpks, locstrainpeak] = findpeaks(xx2,'MinPeakDistance',300,'Threshold',threshold_parameter);

if j == 3
    strainpks = [xx2(1); strainpks];
    
    locstrainpeak = [1; locstrainpeak];
end

%Find local maxima for force
[forcepks , locforcepeak ] = findpeaks(yy2,'MinPeakDistance',300,'Threshold',threshold_parameter);

if plotting == 1
    figure
    plot(locstrainpeak,strainpks,'r*')
    hold on
    plot(xx2)
    xlabel('Strain')
    title('Strain Peaks')
end

if plotting == 1
    figure
    plot(locforcepeak,forcepks,'r*')
    hold on 
    plot(yy2)
    xlabel('Force')
    title('Force Peaks')
end

% Grab just the minimum values for the strain
minimum_values_strain = [];
minimum_values_strain_indices = [];
for i=1:1:length(strainpks)
    if i > 1 && i<length(strainpks)
        if strainpks(i) < strainpks(i-1) && strainpks(i) < strainpks(i+1)
            minimum_values_strain = [minimum_values_strain strainpks(i)];
            minimum_values_strain_indices = [minimum_values_strain_indices locstrainpeak(i)];
        end
    elseif i == length(strainpks)
        if strainpks(i) < strainpks(i-1)
            minimum_values_strain = [minimum_values_strain strainpks(i)];
            minimum_values_strain_indices = [minimum_values_strain_indices locstrainpeak(i)];        
        end    
    end
end


if plotting == 1
    figure
    plot(minimum_values_strain_indices,minimum_values_strain,'r*')
    hold on
    plot(xx2)
    xlabel('Strain')
    title('Minimum Strain Peaks')
end

delta_unloaded = minimum_values_strain;
% Grab just the maximum values for the strain
maximum_values_strain = [];
maximum_values_strain_indices = [];
for i=1:1:length(strainpks)
    if i > 1 && i<length(strainpks)
        if strainpks(i) > strainpks(i-1) && strainpks(i) > strainpks(i+1)
            maximum_values_strain = [maximum_values_strain strainpks(i)];
            maximum_values_strain_indices = [maximum_values_strain_indices locstrainpeak(i)];
        end
    elseif i == length(strainpks)
        if strainpks(i) > strainpks(i-1)
            maximum_values_strain = [maximum_values_strain strainpks(i)];
            maximum_values_strain_indices = [maximum_values_strain_indices locstrainpeak(i)];            
        end 
    elseif i == 1
        if strainpks(i) > strainpks(i+1)
            maximum_values_strain = [maximum_values_strain strainpks(i)];
            maximum_values_strain_indices = [maximum_values_strain_indices locstrainpeak(i)];               
        end       
    end
end
%Save the values of alphaMax that are kept constant during the
%unloading/reloading
alphaMaxUnloading = maximum_values_strain;
if plotting == 1
    figure
    plot(maximum_values_strain_indices,maximum_values_strain,'r*')
    hold on
    plot(xx2)   
    xlabel('Strain')
    title('Maximum Strain Peaks')
end

% Now we want to know we we are in yielding after we have unloaded
% Every pair is the minimum and the next maximum, we want to know
% when we pass the previous maximum for the stretch delimited by the pair

yield_index_array = zeros(1,length(minimum_values_strain_indices)-1);
for i=1:1:length(minimum_values_strain_indices)-1
    yield_index = find(xx2(minimum_values_strain_indices(i):1:maximum_values_strain_indices(i+1))...
                        > maximum_values_strain(i));
    if ~isempty(yield_index)
        yield_index_array(i) = yield_index(1)+minimum_values_strain_indices(i);
    else
        % If yield not found
        yield_index_array(i) = maximum_values_strain_indices(i);
    end

end

if plotting == 1
    figure
    plot(yield_index_array,xx2(yield_index_array),'r*')
    hold on
    plot(xx2)   
    xlabel('Strain')
    title('New yielding points')
end


% Now we can finally get the values for the variable alphaMax, the criteria
% is like this. It is the same to the deformation if we are yielding, it
% doesn't change in the elastic reloading and it is updated again when we
% enter the yielding region

delta = xx2;
alphaMax = xx2;

force = yy2;

delta_plastic_regime = alphaMax;
force_plastic_regime = force;

delta_elastic_regime = cell(length(minimum_values_strain_indices),1);
force_elastic_regime = cell(length(minimum_values_strain_indices),1);
alphaMax_elastic_regime = cell(length(minimum_values_strain_indices),1);

for i=1:1:length(yield_index_array)
    alphaMax(maximum_values_strain_indices(i):1:yield_index_array(i)) = maximum_values_strain(i);
    delta_plastic_regime(maximum_values_strain_indices(i):1:yield_index_array(i)) = 0;
    force_plastic_regime(maximum_values_strain_indices(i):1:yield_index_array(i)) = 0;
    
    if i <= length(minimum_values_strain_indices)
        delta_elastic_regime{i} =xx2(maximum_values_strain_indices(i):1:minimum_values_strain_indices(i));
        alphaMax_elastic_regime{i} = alphaMax(maximum_values_strain_indices(i):1:minimum_values_strain_indices(i));
        force_elastic_regime{i} = force(maximum_values_strain_indices(i):1:minimum_values_strain_indices(i));
    end
end
% Include last unloading curve
delta_plastic_regime(maximum_values_strain_indices(end):1:end) = 0;
force_plastic_regime(maximum_values_strain_indices(end):1:end) = 0;

% Include the alphaMax for the last unloading curve
alphaMax(maximum_values_strain_indices(end):end) = maximum_values_strain(end);

delta_elastic_regime{end} =xx2(maximum_values_strain_indices(end):1:minimum_values_strain_indices(end));
alphaMax_elastic_regime{end} = alphaMax(maximum_values_strain_indices(end):1:minimum_values_strain_indices(end));
force_elastic_regime{end} = force(maximum_values_strain_indices(end):1:minimum_values_strain_indices(end));



%Remove zeros from the delta_plastic_regime and force_plastic_regime. This
%is just removing the unloading/reloading sections.
delta_plastic_regime(delta_plastic_regime == 0) = [];
force_plastic_regime(force_plastic_regime == 0) = [];

% With alphaP I intend do rule out the first section of the yielding curve,
% where there is more curvature.
alphaY = 10*alphaY;
temp = delta_plastic_regime;
% delta_plastic_regime = delta_plastic_regime( temp > alphaY);
% force_plastic_regime = force_plastic_regime ( temp > alphaY );

smooth_parameter = 511;
delta_plastic_regime_filter = smooth(delta_plastic_regime( temp > alphaY),smooth_parameter);
force_plastic_regime_filter = smooth(force_plastic_regime ( temp > alphaY ),smooth_parameter); 

force_plastic_regime_filter = [force_plastic_regime( temp <= alphaY); force_plastic_regime_filter];
delta_plastic_regime_filter = [delta_plastic_regime( temp <= alphaY); delta_plastic_regime_filter];

if plotting == 1
    figure
    plot(delta_plastic_regime_filter,force_plastic_regime_filter)
    hold on
    plot(delta_plastic_regime,force_plastic_regime,'r')
    xlabel('delta')
    ylabel('Force')
    title('Plastic regime')
end

force_plastic_regime = force_plastic_regime_filter;
delta_plastic_regime = delta_plastic_regime_filter;

% Plotting the elastic regime
if plotting == 1
    figure
    for i=1:1:length(delta_elastic_regime)
        plot(delta_elastic_regime{i},force_elastic_regime{i})
        hold on
    end 
    title('Elastic regime')
end

%I'm shifting this array because for each "time step" I want to have the
%current plastic deformation to see where I am. Then, I apply the delta and
%see if I'm yielding or not. If I am, I would update the alphaMax (although
%in this case it just takes the next value in the array
alphaMax = [0; alphaMax(1:end-1)];

if plotting == 1
    figure
    plot(alphaMax)   
    hold on
    plot(xx2,'r-')
    xlabel('Plastic deformation')
    title('AlphaMax')
end

%I'm going to grab only the data from the unloading and not the elastic
%reloading, they seem to be pretty close, but I don't want to make things
%complicated for the curve fitting. I'm going to discard the data from the 
%minimum_values_strain till the yield index.

new_delta = [];
new_force = [];
new_alphaMax = [];

n_minimum_strain = length(minimum_values_strain_indices);
yield_index_array_copy = [1 yield_index_array];
if plotting == 1
    figure
end
for i=1:1:n_minimum_strain
    delta_stretch = delta(yield_index_array_copy(i):1:minimum_values_strain_indices(i));
    force_stretch = force(yield_index_array_copy(i):1:minimum_values_strain_indices(i));
    alphaMax_stretch = alphaMax(yield_index_array_copy(i):1:minimum_values_strain_indices(i));
    new_delta = [new_delta;delta_stretch ];
    new_force = [new_force; force_stretch];
    new_alphaMax = [new_alphaMax; alphaMax_stretch];
    if plotting == 1
        plot(delta_stretch,force_stretch) 
        hold on
    end
end
if plotting == 1
    xlabel('Strain')
    ylabel('Force')
    title('Cropped contact law')
end


end