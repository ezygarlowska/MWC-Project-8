function st = stdevOrbital(T,Hrms,k,h)
%STDEVORBITAL: compute standard deviation orbital motion
%   ST = STDEVORBITAL(T, HRMS, K, H) computes the standard deviation
%   of the near-bed orbital velocity using the wave period T, the wave
%   height HRMS, the wave number K, and the total water depth H (which
%   includes wave set-up and tidal influences).

st = (2*pi/T)*Hrms/(2*sqrt(2)*sinh(k*h));
