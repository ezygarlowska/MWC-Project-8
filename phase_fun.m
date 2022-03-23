%this function computes phase velocity from given period T and water depth h
function [c] = phase_fun(T,h)
     k = k_fun(T,h);
     L = 2*pi/k;
     c = L/T;
end