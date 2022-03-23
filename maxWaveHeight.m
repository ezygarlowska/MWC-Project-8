function Hmax = maxWaveHeight(k,gamma,h)
%MAXWAVEHEIGHT compute maximum wave height
%   HMAX = MAXWAVEHEIGHT(K, GAMMA, H) computes the maximum wave height
%   which is used in Battjes-Jansen. MAXWAVEHEIGHT uses as input
%   parameters the wave number K, the breaker parameter GAMMA and the
%   weighted or total water depth H.

Hmax = (0.88/k)*tanh(gamma*k*h/0.88);