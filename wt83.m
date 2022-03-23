function uv = wt83(st,v)
% uv = wt83(st,v)
% this function computes the wt83 parameterization for the velocity moment
% Input
%   ST   standard deviation of the orbital velocity
%   V    alongshore velocity
% Output
%   UV   velocity moemnt


uv = st.*v.*sqrt(1.16^2 + (v./st).^2);
