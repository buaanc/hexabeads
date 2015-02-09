function [delta, alphaMax, force,new_delta,new_force, new_alphaMax] = reading_data(plotting)

M = csvread('Loading-Unloading_to_6kN.csv',1,7,[1 7 11870 8]);    

if plotting == 1
    plot(M(:,1),M(:,2))
    hold on
end
smooth_parameter = 50;
yy2 = smooth(M(:,2),smooth_parameter);
xx2 = smooth(M(:,1),smooth_parameter); 
if plotting == 1
    plot(xx2,yy2, 'r')

    xlabel('Strain')
    ylabel('Force')
    legend('Experimental data','Smoothed experimental data')
end

threshold_parameter = 0;
% Find local minima for strain
[strainpks, locstrainpeak] = findpeaks(xx2,'MinPeakDistance',300,'Threshold',threshold_parameter);

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
    end
end
%Add the last minimum
minimum_values_strain = [minimum_values_strain strainpks(end)];
minimum_values_strain_indices = [minimum_values_strain_indices locstrainpeak(end)];


if plotting == 1
    figure
    plot(minimum_values_strain_indices,minimum_values_strain,'r*')
    hold on
    plot(xx2)
    xlabel('Strain')
    title('Minimum Strain Peaks')
end
% Grab just the maximum values for the strain
maximum_values_strain = [];
maximum_values_strain_indices = [];
for i=1:1:length(strainpks)
    if i > 1 && i<length(strainpks)
        if strainpks(i) > strainpks(i-1) && strainpks(i) > strainpks(i+1)
            maximum_values_strain = [maximum_values_strain strainpks(i)];
            maximum_values_strain_indices = [maximum_values_strain_indices locstrainpeak(i)];
        end
    end
end
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
    yield_index_array(i) = yield_index(1)+minimum_values_strain_indices(i);

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

for i=1:1:length(yield_index_array)
    alphaMax(maximum_values_strain_indices(i):1:yield_index_array(i)) = maximum_values_strain(i);
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
figure
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
xlabel('Strain')
ylabel('Force')
title('Cropped contact law')


end