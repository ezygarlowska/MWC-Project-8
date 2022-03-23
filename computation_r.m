function [r] = computation_r(S,A)
% computation of the non-linearity parameter r (see Eq. 4, Ruessink et al. 2012)
% input   S skewness
%         A asymmetry
% output
%         r non-linearity parameter

% total non-linearity B (Eq. 7)
B = sqrt(S.^2+A.^2);   
% computation of r (Eq. 11)
r = B2r(B);
return
