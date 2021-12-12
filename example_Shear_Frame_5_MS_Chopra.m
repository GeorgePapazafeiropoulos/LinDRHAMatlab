%% Five-storey shear frame dynamic analysis with modal superposition (Chopra, 2019)

%% Statement of the problem
% * *Chopra (2019), Section 13.2.6:* Consider the five-story shear frame of
% Fig. 12.8.1, subjected to the El Centro ground motion. 
% * *Chopra (2019), Section 13.2.6:* The structure is subjected to the El
% Centro ground motion (Chopra (2012), Fig. 6.1.4). The lumped mass $$m_j =
% m = 45 Mg$ (=0.45kN-sec^2/cm) at each floor, the lateral stiffness of each story is
% $$k_j = k = 54.82 kN/cm.$, and the height of each story is 4 m. The
% damping ratio for all natural modes is $$\mathrm{\zeta_n} = 0.05$.
%
% <<Chopra1281.png>>
%
%% Initialization of structural input data
% Set the storey height of the structure in m.
h=4;
%%
% Set the number of degrees of freedom of the structure, which is equal to
% the number of its storeys.
nDOFs=5;
%%
% Set the lateral stiffness of each storey in N/m.
k=5.482e6;
%%
% Set the lumped mass at each floor in kg.
m=45e3;
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in N/m.
K=k*(diag([2*ones(nDOFs-1,1);1])+diag(-ones(nDOFs-1,1),1)+diag(-ones(nDOFs-1,1),-1));
%%
% Calculate the mass matrix of the structure.
M=m*eye(nDOFs);
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
%% Dynamic Response History Analysis (DRHA) with modal superposition
% Define time
t=dt*(0:(numel(xgtt)-1));
%%
% Perform DRHA analysis for each separate eigenmode to calculate the
% contribution of each eigenmode to the various frame responses

% Initialize
FBeig=cell(5,1);
Feig=cell(5,1);
Ueig=cell(5,1);
MBeig=cell(5,1);
for i=1:5
    % DRHA analysis
    [U,~,~,f] = LDRHA_MS_MDOF(K,M,r,dt,xgtt,ksi,AlgID,u0,ut0,rinf,i);
    % Store the 5th storey shear time history (5th DOF) for each eigenmode
    FBeig{i}=sum(f,1);
    % Store the 5th storey shear time history (5th DOF) for each eigenmode
    Feig{i}=f(5,:);
    % Store the roof displacement time history (5th DOF) for each eigenmode
    Ueig{i}=U(5,:);
    % Base moment time history for each eigenmode
    MBeig{i}=sum(f.*repmat((h:h:5*h)',1,size(f,2)),1);
end
%%
% Perform DRHA analysis for all eigenmodes to calculate the contribution of
% all eigenmodes to the various frame responses

% Eigenmodes that are superposed
eigInd=[1;2;3;4;5];
% DRHA analysis
[U,~,~,f] = LDRHA_MS_MDOF(K,M,r,dt,xgtt,ksi,AlgID,u0,ut0,rinf,eigInd);
% Base shear time history for all eigenmodes
FBeig{6}=sum(f,1);
% 5th storey shear time history (5th DOF) for all eigenmodes
Feig{6}=f(5,:);
% Roof displacement time history (5th DOF) for all eigenmodes
Ueig{6}=U(5,:);
% Base moment time history for all eigenmodes
MBeig{6}=sum(f.*repmat((h:h:5*h)',1,size(f,2)),1);

%% Roof displacement time history
% Plot the contribution of each eigenmode to the roof displacement time
% history. Convert displacements from m to cm.
FigHandle=figure('Name','Roof displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 750]);
for i=1:5
    subplot(nDOFs,1,i)
    plot(t,100*Ueig{i},'LineWidth',1.,'Marker','.',...
        'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
    grid on
    xlim([0,15])
    ylim([-20,20])
    ylabel('U5 (cm)','FontSize',10);
    title(['Mode ',num2str(i),', maxU5=',num2str(max(abs(100*Ueig{i}))),' cm'],...
        'FontSize',10)
end
xlabel('Time (sec)','FontSize',10);
%%
% Plot the contribution of all eigenmodes to the roof displacement time
% history. Convert displacements from m to cm.
FigHandle=figure('Name','Roof displacement','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,100*Ueig{6},'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,15])
ylim([-20,20])
xlabel('Time (sec)','FontSize',10);
ylabel('U5 (cm)','FontSize',10);
title(['All modes, maxU5=',num2str(max(abs(100*Ueig{6}))),' cm'],...
    'FontSize',10)
%%
% Verify with Figure 13.2.8 (left) of Chopra (2019)
%
% <<Chopra1328L.png>>
%
%% Fifth-story shear time history
% Plot the contribution of each eigenmode to the fifth-story shear time
% history. Convert forces from N to kN.
FigHandle=figure('Name','Fifth-story shear','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 750]);
for i=1:5
    subplot(nDOFs,1,i)
    plot(t,Feig{i}/1e3,'LineWidth',1.,'Marker','.',...
        'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
    grid on
    xlim([0,15])
    ylim([-175,175])
    ylabel('V5 (kN)','FontSize',10);
    title(['Mode ',num2str(i),', maxV5=',num2str(max(abs(Feig{i}/1e3))),' kN'],...
        'FontSize',10)
end
xlabel('Time (sec)','FontSize',10);
%%
% Plot the contribution of all eigenmodes to the fifth-story shear time
% history. Convert forces from N to kN.
FigHandle=figure('Name','Fifth-story shear','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,Feig{6}/1e3,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,15])
ylim([-175,175])
xlabel('Time (sec)','FontSize',10);
ylabel('V5 (kN)','FontSize',10);
title(['All modes, maxV5=',num2str(max(abs(Feig{6}/1e3))),' kN'],...
    'FontSize',10)
%%
% Verify with Figure 13.2.7 (right) of Chopra (2019)
%
% <<Chopra1327R.png>>
%
%% Base shear time history
% Plot the contribution of each eigenmode to the base shear time history.
% Convert forces from N to kN.
FigHandle=figure('Name','Base shear','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 750]);
for i=1:5
    subplot(nDOFs,1,i)
    plot(t,FBeig{i}/1e3,'LineWidth',1.,'Marker','.',...
        'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
    grid on
    xlim([0,15])
    ylim([-350,350])
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
plot(t,FBeig{6}/1e3,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,15])
ylim([-350,350])
xlabel('Time (sec)','FontSize',10);
ylabel('Vb (kN)','FontSize',10);
title(['All modes, maxVb=',num2str(max(abs(FBeig{6}/1e3))),' kN'],...
    'FontSize',10)
%%
% Verify with Figure 13.2.7 (left) of Chopra (2019)
%
% <<Chopra1327L.png>>
%
%% Base moment time history
% Plot the contribution of each eigenmode to the base moment time history.
% Convert moments from Nm to kNm.
FigHandle=figure('Name','Base moment','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 750]);
for i=1:5
    subplot(nDOFs,1,i)
    plot(t,MBeig{i}/1e3,'LineWidth',1.,'Marker','.',...
        'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
    grid on
    xlim([0,15])
    ylim([-4000,4000])
    ylabel('Mb (kNm)','FontSize',10);
    title(['Mode ',num2str(i),', maxMb=',num2str(max(abs(MBeig{i}/1e3))),' kNm'],...
        'FontSize',10)
end
xlabel('Time (sec)','FontSize',10);
%%
% Plot the contribution of all eigenmodes to the base moment time history.
% Convert moments from Nm to kNm.
FigHandle=figure('Name','Base moment','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 200]);
plot(t,MBeig{6}/1e3,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
grid on
xlim([0,15])
ylim([-4000,4000])
xlabel('Time (sec)','FontSize',10);
ylabel('Mb (kNm)','FontSize',10);
title(['All modes, maxMb=',num2str(max(abs(MBeig{6}/1e3))),' kNm'],...
    'FontSize',10)
%%
% Verify with Figure 13.2.8 (right) of Chopra (2019)
%
% <<Chopra1328R.png>>
%

%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
