function C = CDM(K,M,ksi)
%
% Classical Damping Matrix of a MDOF system
%
% C = CDM(K,M,ksi)
%
% Description
%     Construct the classical damping matrix of a MDOF system by
%     superposition of modal damping matrices, based on the ratio of
%     critical damping, KSI, given its stiffness and mass matrices, K and M
%     respectively. If the damping matrix of a linear system satisfies the
%     identity c*m^(?1)*k=k*m^(?1)*c the system is said to possess
%     classical damping. This function contains the preferred model for
%     nonlinear RHA of buildings because it eliminates the spurious damping
%     forces that may be encountered in the Rayleigh or Caughey damping
%     models.
%
% Input parameters
%     K [double(:ndofs x :ndofs)] is the stiffness matrix of the system
%         containing only the free degrees of freedom.
%     M [double(:ndofs x :ndofs)] is the mass matrix of the structure
%         containing only the free degrees of freedom.
%     KSI [double(:ndofs x 1)]: ratio of critical damping of the MDOF
%         system for each eigenmode
%
% Output parameters
%     C [double(:ndofs x :ndofs)] is the damping matrix of the structure
%         containing only the free degrees of freedom.
%
%__________________________________________________________________________
% Copyright (c) 2015-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

ndofs=size(K,1);
[Eigvec,Eigval] = eig(K,M);
C=zeros(ndofs);
for i=1:ndofs
    c=2*ksi(i)*sqrt(Eigval(i,i))*(((M*Eigvec(:,i))*Eigvec(:,i)')*M);
    C=C+c;
end
end

