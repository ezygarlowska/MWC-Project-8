% Computation of the sediment transport using the SANTOSS model for a
% coarse sand and skewed wave

% -------------------------------------------------
%                   INITIALISATION
% -------------------------------------------------
clear all
close all
clc


%%% Input parameters for the SANTOSS model

% Wave characteristics
Uw = 1.2;   % orbital velocity amplitude   [m/s]
T = 7;      % period [s]

% Wave shape parameters 
r = 0:0.15:0.6;  % define "r" vector from 0 to 0.6 with step size of 0.15
Nr = length(r);  % number of elements of r
PHI = -pi/2;     % skewed wave 

% Sediment characteristics
D50 = 0.3;  % D50 in mm
D90 = D50;  % D90 in mm
Rhos = 2650; % sediment density in kg/m^3

%%% Initialisation of the output vectors
Qsx = zeros(1,Nr);  % net sediment in the x-direction (m2/s)
Qsy = zeros(1,Nr);  % net sediment in the y-direction (m2/s)
Occ = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the crest
Oct = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otc = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ott = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the trough

% -------------------------------------------------
%          COMPUTATION OF THE SED TRANSPORT
% -------------------------------------------------

for rI = 1:Nr     % loop on the different values of r considered
        
        % 1- computation of the time-series of orbital velocity 
        % ???
        
        % 2- computation of the velocity skewness R and the acceleration skewness beta
        % ???
        
        % 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
        % ???
        
        % 4- sediment transport calculation
        [Qsx(rI) Qsy(rI) Occ(rI) Oct(rI) Ott(rI) Otc(rI)] = SANTOSSmodel(D50,D90,Rhos,T,Urms,R,beta,0,0);
end;


% -------------------------------------------------
%                 VISUALISATION
% -------------------------------------------------
% ???
% ???
