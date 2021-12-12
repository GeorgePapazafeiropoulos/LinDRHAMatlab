%% Two-storey shear frame dynamic analysis with modal superposition (NPTEL)

%% Statement of the problem
% * This example comes from the Introduction to Earthquake Engineering - Web
% course of NPTEL (National Programme on Technology Enhanced Learning),
% Chapter 3, Dynamics of Earthquake Analysis, Example 3.6.
% * A two-story building is modeled as 2-DOF system and rigid floors as
% shown in the following figure. Determine the top floor maximum
% displacement and base shear due to El-Centro, 1940 earthquake ground
% motion. Take the inter-story stiffness, k =197.392 × 103 N/m, the floor
% mass, m = 2500 kg and damping ratio as 2%.
%
% <<NPTEL312.png>>
%
%% Initialization of structural input data
% Set the number of degrees of freedom of the structure, which is equal to
% the number of its storeys.
nDOFs=2;
%%
% Set the lateral stiffness of each storey in N/m.
k=197.392e3*[2;1];
%%
% Set the lumped mass at each floor in kg.
m=2500*[2;1];
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in N/m.
K=diag([0;k])+diag([k;0])+diag(-k,1)+diag(-k,-1);
K=K(2:end,2:end);
%%
% Calculate the mass matrix of the structure.
M=diag(m);
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
r=ones(nDOFs,1);
%% Load earthquake data
% Earthquake acceleration time history of the El Centro earthquake (El
% Centro, 1940, El Centro Terminal Substation Building)
D=load('elcentro.dat');
dt=D(2,1)-D(1,1);
xgtt=9.81*D(:,2);
%%
% Set the critical damping ratio
% ($$\mathrm{\xi}=0.05$)
ksi=0.02;
%%
% Time integration algorithm
AlgID='U0-V0-Opt';
%%
% Initial displacement
u0=zeros(nDOFs,1);
%%
% Initial velocity
ut0=zeros(nDOFs,1);
%%
% Minimum absolute value of the eigenvalues of the amplification matrix
rinf=1;
%% Dynamic Response History Analysis (DRHA) with modal superposition 
% Define time
t=dt*(0:(numel(xgtt)-1));
%%
% Perform DRHA analysis for each separate eigenmode to calculate the
% contribution of each eigenmode to the various frame responses

% Initialize
FBeig=cell(nDOFs+1,1);
Ueig=cell(nDOFs+1,1);
for i=1:nDOFs
    % DRHA analysis
    [U,~,~,f] = LDRHA_MS_MDOF(K,M,r,dt,xgtt,ksi,AlgID,u0,ut0,rinf,i);
    % Store the 2nd storey shear time history (2nd DOF) for each eigenmode
    FBeig{i}=sum(f,1);
    % Store the roof displacement time history (2nd DOF) for each eigenmode
    Ueig{i}=U(2,:);
end
%%
% Perform DRHA analysis for all eigenmodes to calculate the contribution of
% all eigenmodes to the various frame responses

% Eigenmodes that are superposed
eigInd=(1:nDOFs)';
% DRHA analysis
[U,~,~,f] = LDRHA_MS_MDOF(K,M,r,dt,xgtt,ksi,AlgID,u0,ut0,rinf,eigInd);
% Base shear time history for all eigenmodes
FBeig{nDOFs+1}=sum(f,1);
% Roof displacement time history (2nd DOF) for all eigenmodes
Ueig{nDOFs+1}=U(2,:);

%% Roof displacement time history
% Plot the contribution of each eigenmode to the roof displacement time
% history.
FigHandle=figure('Name','Roof displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 750]);
for i=1:nDOFs
    subplot(nDOFs,1,i)
    plot(t,Ueig{i},'LineWidth',1.,'Marker','.',...
        'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
    grid on
    xlim([0,30])
    ylim([-0.25,0.25])
    ylabel('U2 (m)','FontSize',10);
    title(['Mode ',num2str(i),', maxU2=',num2str(max(abs(Ueig{i}))),' m'],...
        'FontSize',10)
end
xlabel('Time (sec)','FontSize',10);
%%
% Plot the contribution of all eigenmodes to the roof displacement time
% history. Convert displacements from m to cm.
FigHandle=figure('Name','Roof displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,Ueig{nDOFs+1},'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-0.25,0.25])
xlabel('Time (sec)','FontSize',10);
ylabel('U2 (m)','FontSize',10);
title(['All modes, maxU2=',num2str(max(abs(Ueig{nDOFs+1}))),' m'],...
    'FontSize',10)
%%
% Verify with Figure 3.13 of Chapter 3, Dynamics of Earthquake Analysis,
% Example 3.6 of NPTEL
%
% <<NPTEL313.png>>
%
%% Base shear time history
% Plot the contribution of each eigenmode to the base shear time history.
% Convert forces from N to kN.
FigHandle=figure('Name','Base shear','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 750]);
for i=1:nDOFs
    subplot(nDOFs,1,i)
    plot(t,FBeig{i}/1e3,'LineWidth',1.,'Marker','.',...
        'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
    grid on
    xlim([0,30])
    ylim([-50,50])
    ylabel('Vb (kN)','FontSize',10);
    title(['Mode ',num2str(i),', maxVb=',num2str(max(abs(FBeig{i}/1e3))),' kN'],...
        'FontSize',10)
end
xlabel('Time (sec)','FontSize',10);
%%
% Plot the contribution of all eigenmodes to the base shear time history.
% Convert forces from N to kN. 
FigHandle=figure('Name','Base shear','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,FBeig{nDOFs+1}/1e3,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-50,50])
xlabel('Time (sec)','FontSize',10);
ylabel('Vb (kN)','FontSize',10);
title(['All modes, maxVb=',num2str(max(abs(FBeig{nDOFs+1}/1e3))),' kN'],...
    'FontSize',10)
%%
% Verify with Figure 3.14 of Chapter 3, Dynamics of Earthquake Analysis,
% Example 3.6 of NPTEL
%
% <<NPTEL314.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
