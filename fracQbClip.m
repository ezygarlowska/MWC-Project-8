function Qb = fracQbClip(Hrms,Hmax,Qb1)
%FRACQBCLIP compute fraction of breaking waves
%   QB = FRACQBCLIP(HRMS, HMAX, QB1) iteratively computes a solution
%   for the fraction of breaking waves QB, using the wave height HRMS,
%   the maximum wave height HMAX, and a first guess for QB, QB1. In
%   this fashion the fraction is the clipped part of the Rayleigh
%   distribution.

tol = 1e-3;
maxiter = 15;
waveratio = -(Hrms/Hmax)^2;

if (Hrms >= Hmax) || isnan(Hrms) || isnan(Hmax)
   Qb = -waveratio;
else
   for i = 1:maxiter
      Qb2 = Qb1;
      fQb = Qb2 - exp((1-Qb2)/waveratio);
      dfQbdQb = 1 + exp((1-Qb2)/waveratio)/waveratio;
      Qb1 = Qb2 - fQb/dfQbdQb;
      if Qb1 < 0, Qb1 = 0; end
      if Qb1 > 1,  Qb1 = 0.95; end
      if abs(Qb2-Qb1) < tol, Qb = Qb1; break, end
   end
end