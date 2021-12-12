%% Classical damping matrix with superposition (Chopra, 2019)

%% Statement of the problem
% * *Chopra (2019), Example 11.1:* The properties of a three-story shear
% building are given in Fig. E11.1.
% * *Chopra (2019), Example 11.3:* Determine a damping matrix for the
% system of Fig. E11.1 by superposing the damping matrices for the first
% two modes, each with $$\mathrm{\zeta_n} = 0.05$.
%
% <<ChopraE111.png>>
%
%% Initialization of structural input data
% Set the lateral stiffness of each storey in N/m.
k=[1e8;1e8;1e8];
%%
% Set the lumped mass at each floor in kg.
m=[2e5;2e5;1e5];
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in N/m.
K=diag([k;0])+diag([0;k])-diag(k,1)-diag(k,-1);
K=K(2:end,2:end);
%%
% Calculate the mass matrix of the structure.
M=diag(m);
%%
% Set the damping ratios for the various eigenmodes
ksi=[0.05;0.05;0];
%% Construct the damping matrix
% Classical damping matrix with modal superposition.
C = CDM(K,M,ksi);
%%
% Convert to kN-sec/cm.
C=1e-5*C

%%
% Verify with Example 11.3 of Chopra (2019)
%
% <<ChopraE113.png>>
%

%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
