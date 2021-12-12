%% Industrial building dynamic analysis with modal superposition (NPTEL)

%% Statement of the problem
% * This example comes from the Introduction to Earthquake Engineering -
% Web course of NPTEL (National Programme on Technology Enhanced Learning),
% Chapter 3, Dynamics of Earthquake Analysis, Example 3.7.
% * An industrial structure is modeled as 2-DOF system as shown in the
% Figure below. Determine the horizontal and vertical displacement of the
% free end of the structure due to El-Centro, 1940 earthquake ground
% motion. Take EI =80 × 103 N.m2, L= 2m, m1= 100kg and m2= 200kg. The
% damping shall be considered as 2 percent.
%
% <<NPTEL315.png>>
%
%% Initialization of structural input data
% Set the number of degrees of freedom of the structure, which is equal to
% the number of its storeys.
nDOFs=2;
%%
% Flexural rigidity
EI=80e3;
%%
% Length
L=2;
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in N/m.
K=6*EI/(7*L^3)*[8,-3;-3,2];
%%
% Calculate the mass matrix of the structure.
M=[300,0;0,200];
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
r=[1;0];
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
UVeig=cell(nDOFs+1,1);
UHeig=cell(nDOFs+1,1);
for i=1:nDOFs
    % DRHA analysis
    [U,~,~,~] = LDRHA_MS_MDOF(K,M,r,dt,xgtt,ksi,AlgID,u0,ut0,rinf,i);
    % Store the vertical displacement for each eigenmode
    UVeig{i}=U(2,:);
    % Store the horizontal displacement for each eigenmode
    UHeig{i}=U(1,:);
end
%%
% Perform DRHA analysis for all eigenmodes to calculate the contribution of
% all eigenmodes to the various frame responses

% Eigenmodes that are superposed
eigInd=(1:nDOFs)';
% DRHA analysis
[U,~,~,f] = LDRHA_MS_MDOF(K,M,r,dt,xgtt,ksi,AlgID,u0,ut0,rinf,eigInd);
% Vertical displacement for all eigenmodes
UVeig{nDOFs+1}=U(2,:);
% Horizontal displacement for all eigenmodes
UHeig{nDOFs+1}=U(1,:);

%% Horizontal displacement time history
% Plot the contribution of each eigenmode to the horizontal displacement
% time history.
FigHandle=figure('Name','Horizontal displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 750]);
for i=1:nDOFs
    subplot(nDOFs,1,i)
    plot(t,UHeig{i},'LineWidth',1.,'Marker','.',...
        'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
    grid on
    xlim([0,30])
    ylim([-0.05,0.05])
    ylabel('U (m)','FontSize',10);
    title(['Mode ',num2str(i),', maxU=',num2str(max(abs(UHeig{i}))),' m'],...
        'FontSize',10)
end
xlabel('Time (sec)','FontSize',10);
%%
% Plot the contribution of all eigenmodes to the horizontal displacement
% time history.
FigHandle=figure('Name','Horizontal displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,UHeig{nDOFs+1},'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-0.05,0.05])
xlabel('Time (sec)','FontSize',10);
ylabel('U (m)','FontSize',10);
title(['All modes, maxU=',num2str(max(abs(UHeig{nDOFs+1}))),' m'],...
    'FontSize',10)
%%
% Verify with Figure 3.16 of Chapter 3, Dynamics of Earthquake Analysis,
% Example 3.7 of NPTEL
%
% <<NPTEL316.png>>
%
%% Vertical displacement time history
% Plot the contribution of each eigenmode to the vertical displacement time
% history.
FigHandle=figure('Name','Vertical displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 750]);
for i=1:nDOFs
    subplot(nDOFs,1,i)
    plot(t,UVeig{i},'LineWidth',1.,'Marker','.',...
        'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
    grid on
    xlim([0,30])
    ylim([-0.07,0.07])
    ylabel('V (m)','FontSize',10);
    title(['Mode ',num2str(i),', maxV=',num2str(max(abs(UVeig{i}))),' m'],...
        'FontSize',10)
end
xlabel('Time (sec)','FontSize',10);
%%
% Plot the contribution of all eigenmodes to the vertical displacement time
% history.
FigHandle=figure('Name','Vertical displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,UVeig{nDOFs+1},'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,30])
ylim([-0.07,0.07])
xlabel('Time (sec)','FontSize',10);
ylabel('V (m)','FontSize',10);
title(['All modes, maxV=',num2str(max(abs(UVeig{nDOFs+1}))),' m'],...
    'FontSize',10)
%%
% Verify with Figure 3.17 of Chapter 3, Dynamics of Earthquake Analysis,
% Example 3.7 of NPTEL
%
% <<NPTEL317.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
