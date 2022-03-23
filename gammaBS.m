function gamma = gammaBS(Hrms0,k0)
%GAMMABS compute breaking parameter gamma
%   GAMMA = GAMMABS(HRMS0, K0) computes the breaker parameter GAMMA
%   according to Battjes and Stive (1985). GAMMABS uses as input the
%   offshore wave height HRMS0, and the offshore wave number K0.

s0 = Hrms0 / (2*pi/k0);
gamma = 0.5 + 0.4*tanh(33*s0);