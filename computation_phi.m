function [phi] = computation_phi(S,A)
% Computation of the phase phi (see Eq. 4, Ruessink et al 2012)
% input   S skewness
%         A asymmetry
% output
%         phi phase in rad

psi = atan2(A,S);       % Eq. 8, Ruessink et al 2012
phi = -psi - pi/2;      % Eq. 12

end