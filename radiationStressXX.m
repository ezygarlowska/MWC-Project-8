function Sxx = radiationStressXX(n,dir,Ew,Er)
%RADIATIONSTRESSXX compute radiation stress.
%   SXX = RADIATIONSTRESSXX(N, DIR, EW, ER) returns the component of
%   the radiation stress in the cross-shore direction due to the
%   gradient in radiation stress from wave motion, and the
%   wave-averaged roller momentum flux. RADIATIONSTRESSXX uses as
%   input the proportion between the wave group velocity and the
%   celerity N, the direction in which the wave is traveling relative
%   to the shore normal in degrees DIR, the wave energy EW, and the
%   roller energy ER.

Sxx = (n + n*cos(dir)^2 -0.5)*Ew + ... % wave part
      2*Er*cos(dir)^2; % roller part