clear all
clc


N_beads = 3;

dynamic = 1;

alpha = 1.5;
%Big Beads
R  = 0.004763;
E = 115E9;
rho = 8500;
nu = 0.3;

%Wall
Ewall = 3.1E9;
rhowall = 8500;
nuwall = 0.30;
%Striker
Rstriker = 0.004763;
Estriker = 115E9;
nustriker = 0.3;
rhostriker = 8500;
StrikerVelocity = 0.9046;
StrikerVelocityDOF = 1;

finN = fopen('input_Nodes.txt', 'w');
finE = fopen('input_Elements.txt', 'w');
finC = fopen('input_Constraints.txt', 'w');
finI = fopen('input_InitialCondition.txt', 'w'); 
finA = fopen('input_TargetArea.txt', 'w'); 


includeintruders = true;

% =========================================================================
% Big Beads Nodal Information
% =========================================================================

Y = zeros(1,N_beads);
X = 0:2*R:(N_beads - 1)*(2*R);

Evec = E*ones(1,N_beads);
nuvec = nu*ones(1,N_beads);
rhovec = rho*ones(1,N_beads);
Rvec = R*ones(1,N_beads);
type = 1*ones(1,N_beads);
design = 0*ones(1,N_beads);
constrained = 0*ones(1,N_beads);

%Only the last bead as design
design(1:end) = 1;

NumberofNodes = N_beads;
% =========================================================================
% Big Beads Connectivity
% =========================================================================
Connectivity = zeros(N_beads-1,2);
Connectivity(:,1) = 1:1:N_beads-1;
Connectivity(:,2) = 2:1:N_beads;
NumberofElements = N_beads-1;
% =========================================================================
% Adding Wall Constraints
% =========================================================================
% Fixed DOF = 1 and 2 => constflag = 3;
% Fixed DOF = 1       => constflag = 1;
% Fixed DOF = 2       => constflag = 2;
constflag = 3;
constrained(end) = constflag;
constrained(1:end-1) = 2;
type(end) = 2;
NumberOfConstraints = sum(constrained ~= 0) + sum(constrained == 3);
% =========================================================================
% Striker
% =========================================================================

Strikernode = 1;



peak_load = 20  %Ultimo 16820.1754385965
%total_t is calculated with Raj's paper
% T = pi * I / 2* P             I = 1.27e−1 Ns   P = 20 Kn
I = 1.27e-1;
total_t = pi * I / (2*peak_load*1e3)*1e6 ; %microseconds
timestep = 1;
Nsteps = fix(total_t/timestep);
shape = 3; %1 triangular, 2 sine, 3 Raj's input load
type(Strikernode)   = 3;

if dynamic == 1
    fprintf(finI,'%d\n %d\t %d\t %d\t',1, Strikernode,StrikerVelocityDOF,total_t);
    if shape == 1
        StrikerVelocity = linspace(0,peak_load,fix(Nsteps/2)+1);
        StrikerVelocity_prov = flipdim(StrikerVelocity,2);
        StrikerVelocity = [StrikerVelocity StrikerVelocity_prov(2:end)];
        for j=1:Nsteps+1
            fprintf(finI,'%f\t ',StrikerVelocity(j));
        end
    elseif shape == 2
        StrikerVelocity = peak_load/2 + (peak_load/2)*sin(pi/(Nsteps/2)*[0:1:Nsteps]-pi/2);

        for j=1:Nsteps+1
            fprintf(finI,'%f\t ',StrikerVelocity(j));
        end
    elseif shape == 3 % Raj's input load
        StrikerVelocity = peak_load*sin(pi/total_t * linspace(0,total_t, total_t/timestep));

        for j=1:Nsteps
            fprintf(finI,'%f\t ',StrikerVelocity(j));
        end
    end
else
    fprintf(finI,'%d\n %d\t %d\t %d\t',1, Strikernode,StrikerVelocityDOF,peak_load);
end

% =========================================================================
% Printing nodal info and connectivity
% =========================================================================
fprintf(finN,'%d\n',NumberofNodes);
fprintf(finC,'%d\n',NumberOfConstraints);
for i=1:NumberofNodes
    fprintf(finN,'%d\t %6.9f\t %6.9f\t %6.9f\t %d\t %10.3e\t %10.3e\t  %6.2f\t %6.9f\n ', ...
            i, X(i),Y(i),1.0,type(i),Evec(i),rhovec(i),nuvec(i),Rvec(i));
    if (constrained(i) ~= 0)
        if (constrained(i) == 3)
            fprintf(finC,'%d\t %d\t %f\n ', i,1,0.0);
            fprintf(finC,'%d\t %d\t %f\n ', i,2,0.0);
        else
            fprintf(finC,'%d\t %d\t %f\n ', i,constrained(i),0.0);
        end
    end
end
fprintf(finE,'%d\n',NumberofElements);
for i=1:NumberofElements
    fprintf(finE,'%d\t %d\t %d\t %f\n ', i, Connectivity(i,1),Connectivity(i,2),i);
end

fprintf(finA,'%d\n',1);
fprintf(finA,'%d\t ', NumberofElements);
fprintf(finA,'\n');

,
fclose(finN);
fclose(finE);
fclose(finC);
fclose(finI);
