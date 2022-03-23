function Dw = dissBreakingBJ(alpha,Qb,T,Hmax,rho,g)
%DISSBREAKINGBJ compute dissipation due to wave breaking
%   DW = DISSBREAKINGBJ(ALPHA, QB, T, HMAX) computes the dissipation
%   due to wave breaking according to Battjes-Jansen. The function
%   uses as inputs the breaker parameter ALPHA, the fraction of
%   breaking waves QB, the wave period T, the maximum wave height
%   HMAX, the density of water RHO and the gravitational
%   accelaration G.

Dw = (0.25/T)*rho*g*alpha*Qb*Hmax^2;