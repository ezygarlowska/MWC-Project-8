%This function computes the root mean square error between two vectors of
%same length. First input vector contains modelled values and the second one
%the observed values. 
function [erms]=root_mean_square_error(model,observations)
    n=length(observations); % number of observations
    erms=0; %initialisation
    for i=1:n
       erms=erms+((model(i)-observations(i)).^2); 
    end 
    erms=sqrt((1/n)*erms); 