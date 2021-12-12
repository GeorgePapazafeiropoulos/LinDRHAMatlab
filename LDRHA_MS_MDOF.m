function [U,V,A,F] = LDRHA_MS_MDOF(K,M,r,dt,xgtt,ksi,AlgID,u0,ut0,rinf,eigInd)
%
% Linear Dynamic Response History Analysis (DRHA) of a MDOF system with
% eigenmode superposition
%
% [U,V,A,F] = LDRHA_MS_MDOF(K,M,R,DT,XGTT,KSI,ALGID,U0,UT0,RINF,EIGIND)
%
% Description
%     Determine the time history of the linear structural response of a
%     Multi-DOF (MDOF) system by linear superposition of its modal
%     responses.
%
% Input parameters
%     K [double(:ndofs x :ndofs)] is the stiffness matrix of the system
%         containing only the free degrees of freedom.
%     M [double(:ndofs x :ndofs)] is the mass matrix of the structure
%         containing only the free degrees of freedom.
%     R [double(:ndofs x 1)]: influence vector. It determines the spatial
%         distribution of the effective earthquake forces
%     DT [double(1 x 1)]: uniform time step of acceleration time history
%         XGTT
%     XGTT [double(:nstep x 2)] 2-column vector of the acceleration history
%         of the excitation imposed at the base. The first column contains
%         time and the second column contains acceleration. nstep is the
%         number of time steps of the dynamic response.
%     KSI [double(1 x 1)]: ratio of critical damping of the MDOF system
%     ALGID [char(1 x :inf)]: algorithm to be used for the time
%         integration. It can be one of the following strings for superior
%         optimally designed algorithms:
%             'generalized a-method': The generalized a-method (Chung &
%             Hulbert, 1993)
%             'HHT a-method': The Hilber-Hughes-Taylor method (Hilber,
%             Hughes & Taylor, 1977)
%             'WBZ': The Wood–Bossak–Zienkiewicz method (Wood, Bossak &
%             Zienkiewicz, 1980)
%             'U0-V0-Opt': Optimal numerical dissipation and dispersion
%             zero order displacement zero order velocity algorithm
%             'U0-V0-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) zero order displacement zero order
%             velocity algorithm
%             'U0-V0-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) zero order displacement zero order
%             velocity algorithm
%             'U0-V1-Opt': Optimal numerical dissipation and dispersion
%             zero order displacement first order velocity algorithm
%             'U0-V1-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) zero order displacement first order
%             velocity algorithm
%             'U0-V1-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) zero order displacement first order
%             velocity algorithm
%             'U1-V0-Opt': Optimal numerical dissipation and dispersion
%             first order displacement zero order velocity algorithm
%             'U1-V0-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) first order displacement zero order
%             velocity algorithm
%             'U1-V0-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) first order displacement zero order
%             velocity algorithm
%             'Newmark ACA': Newmark Average Constant Acceleration method
%             'Newmark LA': Newmark Linear Acceleration method
%             'Newmark BA': Newmark Backward Acceleration method
%             'Fox-Goodwin': Fox-Goodwin formula
%     U0 [double(:ndofs x 1)]: initial displacement of the MDOF system
%     UT0 [double(:ndofs x 1)]: initial velocity of the MDOF system
%     RINF [double(1 x 1)]: minimum absolute value of the eigenvalues of
%         the amplification matrix. For the amplification matrix see
%         eq.(61) in Zhou & Tamma (2004).
%     EIGIND [double(:inf x 1)] is the eigenmode indicator. Only the
%         eigenmode numbers that are contained in EIGIND are taken into
%         account for the modal superposition.
%
% Output parameters
%     U ([ndofs x nstep]): displacement time history.
%     V ([ndofs x nstep]): velocity time history.
%     A ([ndofs x nstep]): acceleration time history.
%     F ([ndofs x nstep]): equivalent static force time history.
%
%__________________________________________________________________________
% Copyright (c) 2015-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

% Number of time integration steps
nstep=length(xgtt);
% Number of degrees of freedom
ndofs=size(M,1);
% Number of eignmodes
neig=numel(eigInd);
% Calculate eigenvalues, eigenvectors and their number
[Eigvec,Eigval]=eig(K,M);
Eigvec=Eigvec(:,eigInd);
% Take the eigenvalues in column vector and sort in ascending order of
% eigenfrequency:
D1=diag(Eigval,0);
D1=D1(eigInd);
% Assemble damping matrix
C=zeros(ndofs);
for i=1:neig
    c=2*ksi*sqrt(D1(i))*(((M*Eigvec(:,i))*Eigvec(:,i)')*M);
    C=C+c;
end
% Generalized masses Mn for all eigenmodes from eq.(13.1.5) of Chopra
% (2012).
Mn=diag(Eigvec'*M*Eigvec);
% Ln coefficients from eq.(13.1.5) of Chopra (2012).
Ln=Eigvec'*M*r;
% Gamman coefficients from eq.(13.1.5) of Chopra (2012).
Gamman=Ln./Mn;
% Eigenperiods of the building
omega=D1.^0.5;
% Initial displacements
u0mod=Eigvec\u0;
% Normalization
u0mod=u0mod./Mn;
% Initial velocities
ut0mod=Eigvec\ut0;
% Normalization
ut0mod=ut0mod./Mn;
% Displacements, velocities and accelerations of the response of the
% eigenmodes of the structure for the given earthquake
U=zeros(neig,nstep);
V=zeros(neig,nstep);
A=zeros(neig,nstep);
F=zeros(neig,nstep);
for i=1:neig
    [u,ut,utt,~] = LDRHA_SDOF(omega(i)^2,1,dt,xgtt,ksi,AlgID,u0mod(i),ut0mod(i),rinf);
    U=U+Gamman(i)*Eigvec(:,i)*u';
    V=V+Gamman(i)*Eigvec(:,i)*ut';
    A=A+Gamman(i)*Eigvec(:,i)*utt';
    F=F+Gamman(i)*omega(i)^2*(M*Eigvec(:,i))*u';
end
end

