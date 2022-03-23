function r = empiricalFitForR(B)
% computes based on an empirical formula a first estimation of r for a given
% value of total nonlinearity B (used in function B2r.m)

% parameter for fit
rPar = 0.9305;

% estimate r
r = tanh(rPar*B);

% ready
return