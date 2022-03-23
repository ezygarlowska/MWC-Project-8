%Here we create the funciton to calculate the wave number k (from period T and water depth h):
function [k] = k_fun(T,h)
    beta = 2.4908;
    g = 9.81;
    x = (h*(2*pi/T))/(sqrt(g*h));
    y = x^2*(1-exp(-x^(beta)))^(-1/(beta));
    k = y/h;
end

