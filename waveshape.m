function [u,t] = waveshape(r,phi,Uw,T)
% waveshape produces the waveshape according to the analytical
% formulation given in Abreu et al. (2010), Coastal Engineering, 57,
% 656-667. 
% 
% INPUT
%   r, measure of non-linearity
%   phi, phase (rad)
%   Uw, velocity amplitude, equals (umax-umin)/2  (m/s)
%   T, period (s)
% OUTPUT
%   u, time series of oscillatory velocity (m/s)
%   t, time axis (s), such that zero-upcrossing is at t = 0 
%
% v1, Gerben Ruessink, 13 May 2010

% make time axis of 1000 steps for output
dt = T/1000;
t = 0:dt:T-dt;

% make time axis such that u(0) = 0 m/s
omega = 2*pi/T;
deltaT = (1/omega)*asin(r*sin(phi)/(1+sqrt(1-r^2)));   % Eq. (25)
tC = t - deltaT;

% compute orbital velocity series
f = sqrt(1-r^2);                 % this ensures that (umax-umin)/2 equals Uw
u = Uw*f*(sin(omega*tC) + (r*sin(phi))/(1 + sqrt(1-r^2))) ./ ...
    (1 - r*cos(omega*tC+phi));   % Eq. (7)

return;
