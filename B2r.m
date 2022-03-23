function r = B2r(B)
% B2r links the "total" non-linearity B to the non-linearity parameter r in
% the Abreu et al. parameterisation. Initially, an empirical fit was
% designed, r = tanh(0.9305*B), as embedded in empiricalFitForR.m. There
% is, however, an analytical solution, given by Malarkey and Davies
% (submitted), Free-stream velocity descriptions under waves with skewness
% and asymmetry.
%
% v1, Gerben Ruessink, 11-2-2012

% make initial guess for b (use old fit)
rIni = empiricalFitForR(B);
bIni = rIni / (1 + sqrt(1-rIni^2));

% estimate b 
b = fzero(@(b) B2b(b,B),bIni);

% estimate r
r = fzero(@(r) b2r(r,b),rIni);

% ready
return

function f = B2b(b,B)
f = B - 3*b / sqrt(2*(1-b^2));

function f = b2r(r,b)
f = b - r / (1 + sqrt(1 - r^2));
