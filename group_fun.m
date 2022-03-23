%this function computes group celerity from period T and water depth h
function [c_g] = group_fun(T,h)
    k = k_fun(T,h);
    n = n_fun(k,h);
    c = phase_fun(T,h);
    c_g = n*c;
end