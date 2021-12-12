%% Four-story frame with an appendage, dynamic analysis with direct integration (Chopra, 2019)

%% Statement of the problem
% * *Chopra (2019), Section 13.2.7:* Consider a four-story building with a
% light appendage-a penthouse, a small housing for mechanical equipment, an
% advertising billboard, or the like. This example is presented because it
% brings out certain special response features representative of a system
% with two natural frequencies that are close.
% * *Chopra (2019), Section 13.2.7:* The lumped masses at the first four
% floors are $$m_j = m = 45 Mg$ (=0.45kN-sec^2/cm) and the appendage mass
% is $$m_5 = 0.01m$. The lateral stiffness of each of the first four
% stories is $$k_j = k = 39.431 kN/cm.$, the appendage stiffness $$k_5 =
% 0.0012k$. The height of each story and the appendage is 4 m. The damping
% ratio for all natural modes is $$\mathrm{\zeta_n} = 0.05$. The response
% of this system to the El Centro ground motion is determined
%
%% Initialization of structural input data
% Set the storey height of the structure in m.
h=4;
%%
% Set the number of degrees of freedom of the structure, which is equal to
% the number of its storeys.
nDOFs=5;
%%
% Set the lateral stiffnesses of all storeys in N/m.
k=[3.9431e6;3.9431e6;3.9431e6;3.9431e6;0.0012*3.9431e6];
%%
% Set the lumped mass at each floor in kg.
m=[45e3;45e3;45e3;45e3;0.01*45e3];
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in N/m.
K=diag([0;k])+diag([k;0])-diag(k,1)-diag(k,-1);
K=K(2:end,2:end);
%%
% Calculate the mass matrix of the structure.
M=diag(m);
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
r=ones(5,1);
%% Load earthquake data
% Earthquake acceleration time history of the El Centro earthquake (El
% Centro, 1940, El Centro Terminal Substation Building)
D=load('elcentro.dat');
dt=D(2,1)-D(1,1);
xgtt=9.81*D(:,2);
%%
% Set the critical damping ratio
% ($$\mathrm{\xi}=0.05$)
ksi=0.05;
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
%% Dynamic Response History Analysis (DRHA) with direct integration
% Define time
t=dt*(0:(numel(xgtt)-1));
%%
% Calculate the classical damping matrix of the structure
C = CDM(K,M,ksi*ones(nDOFs,1));
%%
% Perform DRHA analysis
[~,~,~,f] = LDRHA_DI_MDOF(K,C,M,r,dt,xgtt,AlgID,u0,ut0,rinf);
%%
% Base shear time history
FBeig=sum(f,1);
%%
% Appendage shear time history
Feig=f(5,:);
%% Appendage shear time history
% Plot the appendage shear time history. 
FigHandle=figure('Name','Appendage shear','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,Feig/1e3,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,15])
ylim([-9,9])
xlabel('Time (sec)','FontSize',10);
ylabel('V5 (kN)','FontSize',10);
title(['All modes, maxV5=',num2str(max(abs(Feig/1e3))),' kN'],...
    'FontSize',10)
%%
% Verify with Figure 13.2.11 (right) of Chopra (2019)
%
% <<Chopra13211Rtot.png>>
%
%% Base shear time history
% Plot the base shear time history.
FigHandle=figure('Name','Base shear','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,FBeig/1e3,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,15])
ylim([-350,350])
xlabel('Time (sec)','FontSize',10);
ylabel('Vb (kN)','FontSize',10);
title(['All modes, maxVb=',num2str(max(abs(FBeig/1e3))),' kN'],...
    'FontSize',10)
%%
% Verify with Figure 13.2.11 (left) of Chopra (2019)
%
% <<Chopra13211Ltot.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
