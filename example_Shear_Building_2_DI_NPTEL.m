%% Two-storey shear frame, dynamic analysis with direct integration (NPTEL)

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
% Calculate the classical damping matrix of the structure
C = CDM(K,M,ksi*ones(nDOFs,1));
%%
% Perform DRHA analysis 
[U,~,~,f] = LDRHA_DI_MDOF(K,C,M,r,dt,xgtt,AlgID,u0,ut0,rinf);
%%
% Base shear time history for all eigenmodes
FBeig=sum(f,1);
%%
% Roof displacement time history (2nd DOF) for all eigenmodes
Ueig=U(2,:);
%% Roof displacement time history
% Plot the contribution of all eigenmodes to the roof displacement time
% history. Convert displacements from m to cm.
FigHandle=figure('Name','Roof displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,Ueig,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-0.25,0.25])
xlabel('Time (sec)','FontSize',10);
ylabel('U2 (m)','FontSize',10);
title(['All modes, maxU2=',num2str(max(abs(Ueig))),' m'],...
    'FontSize',10)
%%
% Verify with Figure 3.13 of Chapter 3, Dynamics of Earthquake Analysis,
% Example 3.6 of NPTEL
%
% <<NPTEL313tot.png>>
%
%% Base shear time history
% Plot the contribution of all eigenmodes to the base shear time history.
% Convert forces from N to kN. 
FigHandle=figure('Name','Base shear','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,FBeig/1e3,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-50,50])
xlabel('Time (sec)','FontSize',10);
ylabel('Vb (kN)','FontSize',10);
title(['All modes, maxVb=',num2str(max(abs(FBeig/1e3))),' kN'],...
    'FontSize',10)
%%
% Verify with Figure 3.14 of Chapter 3, Dynamics of Earthquake Analysis,
% Example 3.6 of NPTEL
%
% <<NPTEL314tot.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
