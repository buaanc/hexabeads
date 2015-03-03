function hexagonalMATOPT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% input_Nodes.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%
% #Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nodenumber   x   y   z   parttype    designflag   R   rho   nu   R %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E = 115E9;
rho = 8500;

%Wall
% Ewall = 3.1E9;
% rhowall = 1500;
% nuwall = 0.35;
Ewall = 115E9;
rhowall = 8500;
nuwall = 0.3;


% Input number of rows and columns in setup
n = 3; %number of particles along one edge

D=(3/4)*2.54/100; % 3/4 inch to m
r=D/2;
rd=r*(sqrt(2)-1);
Dd=D*(sqrt(2)-1);
x=sqrt(3)*r;

npart = 2*n + (n+n-1) + 6*n;
for i=1:n-2
    npart = npart + 2*(n+i);
end

Nodes = zeros(npart,9);

npack = npart-6*n;
nwall = 6*n-1;


striker_type = 1; % 0 for velocity

if striker_type == 1
    %Loading parameters
    peak_load = 10.1616; %KN High
    %peak_load = 2.843 Low
    %peak_load = 10.01 Medium
    total_t = 48; %microseconds
    timestep = 1.0;
    Nsteps = fix(total_t/timestep);
    shape = 2; %1 triangular, 2 sine
else
    peak_load = 1.616;
end

%%% Numbering %%%
Nodes(:,1) = 1:npart;

%%% x %%%
count = 0;
k=0;
for i=1:n
    for j=1:(n+k)
        count=count+1;
        Nodes(count,2) = (i-1)*x;
    end
    k=k+1;
end

k=n-2;
for i=(n+1):(2*n-1)
    for j=1:(n+k)
        count=count+1;
        Nodes(count,2) = (i-1)*x;
    end
    k=k-1;
end


%wall1
for i=1:n
    count = count+1;
    Nodes(count,2) = -r;
end

%wall2
for i=1:n
    count = count+1;
    Nodes(count,2) = (i-1)*x-1/2*r;
end

count


%wall3
for i=n:(2*n-1)
    count = count+1;
    Nodes(count,2) = (i-1)*x+1/2*r;
end

%wall4
for i=1:n
    count = count+1;
    Nodes(count,2) = (2*n-2)*x+r;
end

%wall5
for i=n-1:1:(2*n-2)
    count = count+1;
    Nodes(count,2) = i*x+1/2*r;
end

%wall6
for i=0:1:n-1
    count = count+1;
    Nodes(count,2) = i*x-1/2*r;
end

%%% y %%%
count = 0;
k=0;
for i=1:n
    for j=1:(n+k)
        count=count+1;
        Nodes(count,3) = (j-1)*D-(i-1)*r;
    end
    k=k+1;
end

k=n-2;
for i=1:n-1
    for j=1:(n+k)
        count=count+1;
        Nodes(count,3) = (j-1)*D-(n-2)*r+(i-1)*r;
    end
    k=k-1;
end
%wall1
for i=1:n
    count = count+1;
    Nodes(count,3) = (i-1)*D;
end

%Place striker in first wall
striker = round(count - n/2);
%wall2
for i=1:n
    count = count+1;
    Nodes(count,3) = (n-1)*D+r+(i-1)*r-(1-sqrt(3)/2)*r;
end

%wall3
for i=1:n
    count = count+1;
    Nodes(count,3) =  (n-1)*D+r+(n-1)*r-(i-1)*r-(1-sqrt(3)/2)*r;
end

%wall4
for i=1:n
    count = count+1;
    Nodes(count,3) = (i-1)*D;
end

%%% Target area
target_nodes = (count-n+1):1:count;

%wall5
for i=1:n
    count = count+1;
    Nodes(count,3) = -(n+1)*r+i*r+(1-sqrt(3)/2)*r;
end

%wall6
for i=1:n
    count = count+1;
    Nodes(count,3) = -i*r+(1-sqrt(3)/2)*r;
end

% count

%%% z %%%
Nodes(:,4)=1;

%%% parttype %%%
Nodes(1:npack,5)=1;
Nodes(npack+1:end,5)=2;

%%% E %%%
Nodes(1:npack,6)= E;
Nodes(npack+1:end,6)= Ewall;

%%% rho %%%
Nodes(1:npack,7) = rho;
Nodes(npack+1:end,7) = rhowall;

%%% nu %%%
Nodes(1:npack,8) = 0.30;
Nodes(npack+1:end,8)=0.35;

%%% R %%%
Nodes(1:npack,9) = r;
Nodes(npack+1:end,9) = r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% input_Elements.txt %%%%%%%%%%%%%%%%%%%%%%%
% #Elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elementnumber   Node1   Node2   alpha %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nelement = 6*n + ((n-1)*(3*n-4)+2*n-2) + (2*(n-1)*(3*n-2));

%%%% DESIGN GROUPS %%%%%

% Each element has its own design variable and does not share it
design_groups = 1:1:nelement;
design_variables = 1:1:nelement;


Elements = zeros(nelement,4);

%%% Numbering %%%
Elements(:,1) = 1:nelement;

%%% Connectivity %%% (packing then walls)
%%% packing first %%
count = 0;
for j=1:n
    for i=(((j-1)*(2*n+j-2)/2)+1):((((j-1)*(2*n+j-2)/2)+1)+n+j-3)
        count = count+1;
        Elements(count,2) = i;
        Elements(count,3) = i+1;
    end
end

for j=1:n-1
    for i=(((j-1)*(2*n+j-2)/2)+1):((((j-1)*(2*n+j-2)/2)+1)+n+j-2)
        count = count+1;
        Elements(count,2) = i;
        Elements(count,3) = i+n+j-1;
        count = count+1;
        Elements(count,2) = i;
        Elements(count,3) = i+n+j;
    end
end

for j=1:n-1
    for i=(npack-(((j-1)*(2*n+j-2)/2)+1)):-1:((npack-(((j-1)*(2*n+j-2)/2)+1))-(n+j-3))
        count = count+1;
        Elements(count,2) = i+1;
        Elements(count,3) = i;
    end
end

for j=1:n-1
    for i=(npack-(((j-1)*(2*n+j-2)/2)+1)):-1:((npack-(((j-1)*(2*n+j-2)/2)+1))-(n+j-2))
        count = count+1;
        Elements(count,2) = i+1;
        Elements(count,3) = i+1-(n+j-1);
        count = count+1;
        Elements(count,2) = i+1;
        Elements(count,3) = i+1-(n+j);
    end
end

%%% walls %%
%left
for i=1:n
        count = count+1;
        Elements(count,2) = i;
        Elements(count,3) = i+npack;
end

%topleft
for i=1:n
        count = count+1;
        Elements(count,2) = i*(2*n+i-1)/2;
        Elements(count,3) = i+npack+n;
end

%topright
for i=1:n
        count = count+1;
        Elements(count,2) = n*(3*n-1)/2+(i-1)*(4*n-i-2)/2;
        Elements(count,3) = i+npack+2*n;
end

%right
for i=1:n
        count = count+1;
        Elements(count,2) = npack-n+i;
        Elements(count,3) = i+npack+3*n;
end

%bottomright
for i=1:n
        count = count+1;
        Elements(count,2) = n*(3*n-1)/2+(i-1)*(4*n-i-2)/2-(2*n-2)+i-1;
        Elements(count,3) = i+npack+4*n;
end

%bottomleft
for i=1:n
        count = count+1;
        Elements(count,2) = i*(2*n+i-1)/2-(n-1+i-1);
        Elements(count,3) = i+npack+5*n;
end

%%% Alpha %%%
Elements(design_groups,4) = design_variables;

%%% Striker
Nodes(striker,5) = 3;
Nodes(striker,end) = r;
Nodes(striker,6) = E;
Nodes(striker,7) = rho;
%Nodes = [Nodes; npart+1, Nodes(striker,2) - D, Nodes(striker,3), 1.0, 3, 1.0, 1.930e+011, 8000,0.30, r];
%npart = npart + 1;

%Elements = [Elements; nelement+1, npart, striker, 1.5];
%nelement = nelement + 1;

scatter(Nodes(:,2),Nodes(:,3))
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NodeFile = fopen('input_Nodes.txt','w+');
fprintf(NodeFile,'%g\n',npart);

for i=1:npart
    fprintf(NodeFile,'%g \t %f \t %f \t %f \t %g \t %e \t %e \t %f \t %e \n',Nodes(i,1),Nodes(i,2),Nodes(i,3),Nodes(i,4),Nodes(i,5),Nodes(i,6),Nodes(i,7),Nodes(i,8),Nodes(i,9));
end
  
fclose(NodeFile);

ElementFile = fopen('input_Elements.txt','w+');
fprintf(ElementFile,'%g\n',nelement);

for i=1:nelement
    fprintf(ElementFile,'%g \t %g \t %g \t %g \n',Elements(i,1),Elements(i,2),Elements(i,3),Elements(i,4));
end




fclose(ElementFile);

ConstraintsFile = fopen('input_Constraints.txt', 'w');
constrained = Nodes(Nodes(:,5) == 2,1);
constrained = repmat(constrained',[2,1]);
constrained = constrained(:);
C_dofs = reshape(repmat([1;2],[1,length(constrained)/2]),[length(constrained),1]);

fprintf(ConstraintsFile,'%g \n',length(constrained));  
for i=1:length(constrained)
    fprintf(ConstraintsFile,'%g \t %g \t %g \n',constrained(i),C_dofs(i),0.0);    
end

for i=1:length(constrained)/2
    constrained = Nodes(Nodes(:,5) == 2,1:4);
    text(constrained(i,2),constrained(i,3),num2str(constrained(i,1)))
end
text(Nodes(end,2),Nodes(end,3),'striker')

fclose(ConstraintsFile);

%%%% Initial condition
IcFile = fopen('input_InitialCondition.txt', 'w'); 


if striker_type == 1
    Strikernode = striker;
    StrikerVelocityDOF = 1;
    fprintf(IcFile,'%d\n %d\t %d\t %d\t',1, Strikernode,StrikerVelocityDOF,Nsteps);
    if shape == 1
        StrikerVelocity = linspace(0,peak_load,fix(Nsteps/2)+1);
        StrikerVelocity_prov = flipdim(StrikerVelocity,2);
        StrikerVelocity = [StrikerVelocity StrikerVelocity_prov(2:end)];
        for j=1:Nsteps+1
            fprintf(IcFile,'%f\t ',StrikerVelocity(j));
        end
    else 
        StrikerVelocity = peak_load/2 + (peak_load/2)*sin(pi/(Nsteps/2)*[0:1:Nsteps]-pi/2);

        for j=1:Nsteps+1
            fprintf(IcFile,'%f\t ',StrikerVelocity(j));
        end
    end

else
    Strikernode = striker;
    StrikerVelocityDOF = 1;
    fprintf(IcFile,'%d\n %d\t %d\t %d\t',1, Strikernode,StrikerVelocityDOF,peak_load);
end
fclose(IcFile);

%%%% Target elements
TargetFile = fopen('input_TargetArea.txt', 'w'); 

TargetElements = zeros(length(target_nodes),1);
% TargetElements = zeros(length(nelement),1);
% TargetElements = 1:1:nelement;
% n_target = nelement;
count = 1;
for i=1:nelement
    if any(Elements(i,2) == target_nodes) || any(Elements(i,3) ==  target_nodes);
        TargetElements(count) = Elements(i,1);
        count = count + 1;
    end
end
n_target = length(target_nodes);

fprintf(TargetFile,'%d\n',n_target);
for i=1:n_target
    fprintf(TargetFile,'%d\t',TargetElements(i));
end


for k=1:nelement
    x_coord = [Nodes(Elements(k,2),2) Nodes(Elements(k,3),2)];
    y_coord = [Nodes(Elements(k,2),3) Nodes(Elements(k,3),3)];
    plot(x_coord,y_coord)
    hold on
end

axis equal

end