clear all
clc

%Loading parameters
peak_load = 22.616; %KN High
%peak_load = 2.843 %Low
%peak_load = 10.01 %Medium

%peak_load = 16.8201754385965  %Ultimo 16820.1754385965
total_t = 10; %microseconds
timestep = 1;
Nsteps = fix(total_t/timestep);
shape = 2; %1 triangular, 2 sine

nbeadx = 3;
nbeady = 2;
alpha = 1.5;
%Big Beads
R  = 0.0095;
E = 115E9;
rho = 8500;
nu = 0.3;
%Intruder Beads
Ri = (sqrt(2.0)-1.0)*R;
Ei = 115E9;
rhoi = 8500;
nui = 0.3;
%Wall
Ewall = 3.1E9;
rhowall = 1500;
nuwall = 0.35;
%Striker
Rstriker = 0.0095;
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

%Get target area for optimization
goal = 'focus_symm';

% =========================================================================
% Big Beads Nodal Information
% =========================================================================
count = 0;
for i=1:nbeadx
    for j=1:nbeady
        count = count + 1;
        X(count) = (i-1)*2*R;
        Y(count) = (j-1)*2*R;
        Evec(count) = E;
        nuvec(count) = nu;
        rhovec(count) = rho;
        Rvec(count) = R;
        type(count) = 1;
        design(count) = 0;
        constrained(count) = 0;
    end
end
NumberofNodes = count;
% =========================================================================
% Big Beads Connectivity
% =========================================================================
count = 0;
for i=1:nbeadx
    for j=1:nbeady-1
        count = count + 1;
        Connectivity(count,1) = (i-1)*nbeady + j;
        Connectivity(count,2) = (i-1)*nbeady + j+1;
    end
end
for i=1:nbeadx-1
    for j=1:nbeady
        count = count + 1;
        Connectivity(count,1) = j + (i-1)*nbeady;
        Connectivity(count,2) = j + (i)*nbeady;
    end
end
NumberofElements = count;
% =========================================================================
% Adding Wall Constraints
% =========================================================================
% Fixed DOF = 1 and 2 => constflag = 3;
% Fixed DOF = 1       => constflag = 1;
% Fixed DOF = 2       => constflag = 2;
constflag = 3;
% =========================================================================
% Left Wall
% =========================================================================
count = NumberofNodes;
NumberOfConstraints = 0;
for i=1:nbeady
    count = count + 1;
    X(count) = -R;
    Y(count) = (i-1)*2*R;
    Evec(count) = Ewall;
    nuvec(count) = nuwall;
    rhovec(count) = rhowall;
    Rvec(count) = R;
    type(count) = 2;
    design(count) = 0;
    constrained(count) = constflag;
    if (constflag == 3)
        NumberOfConstraints = NumberOfConstraints + 2;
    else
        NumberOfConstraints = NumberOfConstraints + 1;
    end
end
NumberofNodes = count;
count = NumberofElements;
for i=1:nbeady
    count = count + 1;
    Connectivity(count,1) = i;
    Connectivity(count,2) = nbeadx*nbeady + i;
end
NumberofElements = count;
% =========================================================================
% Right Wall
% =========================================================================
count = NumberofNodes;
for i=1:nbeady
    count = count + 1;
    X(count) = (nbeadx-1)*2*R + R;
    Y(count) = (i-1)*2*R;    
    Evec(count) = Ewall;
    nuvec(count) = nuwall;
    rhovec(count) = rhowall;
    Rvec(count) = R;
    type(count) = 2;
    design(count) = 0;
    constrained(count) = constflag;
    if (constflag == 3)
        NumberOfConstraints = NumberOfConstraints + 2;
    else
        NumberOfConstraints = NumberOfConstraints + 1;
    end
end
NumberofNodes = count;
count = NumberofElements;
for i=1:nbeady
    count = count + 1;
    Connectivity(count,1) = (nbeadx-1)*nbeady + i;
    Connectivity(count,2) = nbeadx*nbeady + nbeady + i;
end
NumberofElements = count;
% =========================================================================
% Bottom Wall
% =========================================================================
count = NumberofNodes;
for i=1:nbeadx
    count = count + 1;
    X(count) = (i-1)*2*R;
    Y(count) = -R;
    Evec(count) = Ewall;
    nuvec(count) = nuwall;
    rhovec(count) = rhowall;
    Rvec(count) = R;
    type(count) = 2;
    design(count) = 0;
    constrained(count) = constflag;
    if (constflag == 3)
        NumberOfConstraints = NumberOfConstraints + 2;
    else
        NumberOfConstraints = NumberOfConstraints + 1;
    end
end
NumberofNodes = count;
count = NumberofElements;
for i=1:nbeadx
    count = count + 1;
    Connectivity(count,1) = nbeady*(i-1)+1;
    Connectivity(count,2) = nbeady*nbeadx + 2*nbeady +i;
end
NumberofElements = count;
% =========================================================================
% Top Wall
% =========================================================================
count = NumberofNodes;
for i=1:nbeadx
    count = count + 1;
    X(count) = (i-1)*2*R;
    Y(count) = (nbeady-1)*2*R + R;
    Evec(count) = Ewall;
    nuvec(count) = nuwall;
    rhovec(count) = rhowall;
    Rvec(count) = R;
    type(count) = 2;
    design(count) = 0;
    constrained(count) = constflag;
    if (constflag == 3)
        NumberOfConstraints = NumberOfConstraints + 2;
    else
        NumberOfConstraints = NumberOfConstraints + 1;
    end
end
NumberofNodes = count;
count = NumberofElements;
for i=1:nbeadx
    count = count + 1;
    Connectivity(count,1) = nbeady*i;
    Connectivity(count,2) = nbeady*nbeadx + 2*nbeady + nbeadx + i;
end
NumberofElements = count;
% =========================================================================
% Intruders/Small beads nodal information
% =========================================================================
temp  = NumberofNodes;
count = NumberofNodes;
for i=1:nbeadx-1
    for j=1:nbeady-1
        count = count + 1;
        X(count) = R + (i-1)*2*R;
        Y(count) = R + (j-1)*2*R;
        Evec(count) = Ei;
        nuvec(count) = nui;
        rhovec(count) = rhoi;
        Rvec(count) = Ri;
        type(count) = 1;
        design(count) = 1;
        constrained(count) = 0;
    end
end
NumberofNodes = count;
% =========================================================================
% Intruders/Small Connectivity
% =========================================================================
count = NumberofElements;
c1 = 0;
for i=1:nbeadx-1
    for j=1:nbeady-1
        c1 = c1 + 1;
        count = count + 1;
        Connectivity(count,1) = (i-1)*nbeady + j;
        Connectivity(count,2) = temp + c1;
        count = count + 1;
        Connectivity(count,1) = i*nbeady + j;
        Connectivity(count,2) = temp + c1;
        count = count + 1;
        Connectivity(count,1) = (i-1)*nbeady + j + 1;
        Connectivity(count,2) = temp + c1;
        count = count + 1;
        Connectivity(count,1) = i*nbeady + j + 1;
        Connectivity(count,2) = temp + c1;
    end
end
NumberofElements = count;
% =========================================================================
% Striker
% =========================================================================
if  mod(nbeady,2) == 1 % nbeady is odd
    bead = ceil(nbeady/2);
    Strikernode = nbeadx*nbeady + bead;
    X(Strikernode)      = -(Rstriker+R);
    Y(Strikernode)      = (bead-1)*2*R;
    Evec(Strikernode)   = Estriker;
    nuvec(Strikernode)  = nustriker;
    rhovec(Strikernode) = rhostriker;
    Rvec(Strikernode)   = Rstriker;
    type(Strikernode)   = 3;
    design(Strikernode) = 0;
    constrained(Strikernode) = 0;
else % nbeady is even
    bead = nbeady/2;
    h = sqrt(Rstriker^2+2*Rstriker*R);
    NumberofNodes = NumberofNodes + 1;  
    Strikernode = NumberofNodes;
    X(Strikernode) = -h;
    Y(Strikernode) = (bead-1)*2*R + R;
    Evec(Strikernode)   = Estriker;
    nuvec(Strikernode)  = nustriker;
    rhovec(Strikernode) = rhostriker;
    Rvec(Strikernode)   = Rstriker;
    type(Strikernode)   = 3;
    design(Strikernode) = 0;
    constrained(Strikernode) = 0;
    NumberofElements = NumberofElements + 1;
    Connectivity(NumberofElements,1) =  bead;
    Connectivity(NumberofElements,2) =  Strikernode;
    NumberofElements = NumberofElements + 1;
    Connectivity(NumberofElements,1) =  bead+1;
    Connectivity(NumberofElements,2) =  Strikernode;    
end

fprintf(finI,'%d\n %d\t %d\t %d\t',1, Strikernode,StrikerVelocityDOF,Nsteps);
if shape == 1
    StrikerVelocity = linspace(0,peak_load,fix(Nsteps/2)+1);
    StrikerVelocity_prov = flipdim(StrikerVelocity,2);
    StrikerVelocity = [StrikerVelocity StrikerVelocity_prov(2:end)];
    for j=1:Nsteps+1
        fprintf(finI,'%f\t ',StrikerVelocity(j));
    end
else 
    StrikerVelocity = peak_load/2 + (peak_load/2)*sin(pi/(Nsteps/2)*[0:1:Nsteps]-pi/2);
    
    for j=1:Nsteps+1
        fprintf(finI,'%f\t ',StrikerVelocity(j));
    end
end
timeLeft= Nsteps+1:1:800;
[Peak,tMax] = max(StrikerVelocity);
plot([0:1:Nsteps timeLeft],[StrikerVelocity zeros(1,length(timeLeft))])
title('Input Load','Fontsize',15,'Interpreter','latex')
ylabel('$Force: KN$','Interpreter','latex','Fontsize',15)
xlabel('$Time: \mu s$', 'Interpreter','latex','Fontsize',15)


x = [tMax tMax + 2*Nsteps  ];
y = [Peak Peak*2/3];
NormX = x./diff(get(gca,'XLim'));
NormY = y./diff(get(gca,'YLim'));

%a = annotation('textarrow', x_to_norm_v2(x(2),x(1)), y_to_norm_v2(y(2),y(1)) , 'String' , {'Peak = ', num2str(peak_load), 'KN'});
% set(a,'Fontsize',15,'Headstyle','cback3','Interpreter','latex')
% grid on
% =========================================================================
% Printing nodal info and connectivity
% =========================================================================
fprintf(finN,'%d\n',NumberofNodes);
if (mod(nbeady,2) == 1)
    fprintf(finC,'%d\n',NumberOfConstraints-1);
else
    fprintf(finC,'%d\n',NumberOfConstraints);
end
for i=1:NumberofNodes
    fprintf(finN,'%d\t %6.9f\t %6.9f\t %6.9f\t %d\t %d\t %10.3e\t %10.3e\t  %6.2f\t %6.9f\n ', ...
            i, X(i),Y(i),1.0,type(i),design(i),Evec(i),rhovec(i),nuvec(i),Rvec(i));
    if (constrained(i) ~= 0)
        if (constrained(i) == 3)
            fprintf(finC,'%d\t %d\t %f\n ', i,1,0.0);
            fprintf(finC,'%d\t %d\t %f\n ', i,2,0.0);
        else
            fprintf(finC,'%d\t %d\t %f\n ', i,constrained(i),0.0);
        end
    end
end

%%% Adding lateral constraints for the striker for the odd case
if (mod(nbeady,2) == 1)
    fprintf(finC,'%d\t %d\t %f\n ', Strikernode,2,0.0);
end



fprintf(finE,'%d\n',NumberofElements);
for i=1:NumberofElements
    fprintf(finE,'%d\t %d\t %d\t %f\n ', i, Connectivity(i,1),Connectivity(i,2),alpha);
end


% =========================================================================
% Get Target Area for optimization problem
% =========================================================================
W=(nbeadx-1)*2*R;
H=(nbeady-1)*2*R;

TargetNodes = [];
countNodes = 0;
if (strcmp(goal,'disperse'))
    for i=1:NumberofNodes
        if (X(i) == W+R & Y(i) >=0 & Y(i)<= H)
           TargetNodes  = [TargetNodes i];
           countNodes = countNodes + 1;
        end
    end
    TargetArea = zeros(countNodes,1);
    for i=1:countNodes
        node = TargetNodes(i);
        for j=1:NumberofElements
            if (node == Connectivity(j,2))
                TargetArea(i,1) = j;
                break;
            end
        end
    end
elseif (strcmp(goal,'focus_asymm'))  
     for i=1:NumberofNodes
        if (Y(i) ==-R & X(i)<=W/2+3*R & X(i)>=W/2-3*R)
           TargetNodes  = [TargetNodes i];
           countNodes = countNodes + 1;
        end
    end
    TargetArea = zeros(countNodes,1);
    for i=1:countNodes
        node = TargetNodes(i);
        for j=1:NumberofElements
            if (node == Connectivity(j,2))
                TargetArea(i,1) = j;
                break;
            end
        end
    end   
elseif (strcmp(goal,'focus_symm'))   
     for i=1:NumberofNodes
        if ((Y(i) ==-R | Y(i) == H+R) & X(i)<=W/2+3*R & X(i)>=W/2-3*R)
           TargetNodes  = [TargetNodes i];
           countNodes = countNodes + 1;
        end
    end
    TargetArea = zeros(countNodes,1);
    for i=1:countNodes
        node = TargetNodes(i);
        for j=1:NumberofElements
            if (node == Connectivity(j,2))
                TargetArea(i,1) = j;
                break;
            end
        end
    end      
end

fprintf(finA,'%d\n',countNodes);
for i=1:countNodes
    fprintf(finA,'%d\t ', TargetArea(i));
end
fprintf(finA,'\n');

fclose(finN);
fclose(finE);
fclose(finC);
fclose(finI);
fclose(finA);
