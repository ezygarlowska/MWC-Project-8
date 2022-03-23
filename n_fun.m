%this function computes the propagation factor n from the wave number k and water depth h
function [n] = n_fun(k,h)
    n = 0.5*(1+2*k*h/(sinh(2*k*h)));
end