clear; clc; close all;

%% Chapter 9 - Sediment transport modelling

%% Preliminary computations- orbital velocity time-series

T=6; %wave period [s]
Uw=1; %orbital velocity amplitude [m/s]

r=[0 0.6 0.6 0.6]; %is a non-linearity measure
phi=[0 0 -pi/2 -pi/4]; %is a phase

for i=1:4
    [u(:,i),t(:,i)]=waveshape(r(i),phi(i),Uw,T);
    [R(i),Beta(i)]=velocity_skewness_asymmetry(u(:,i),t(:,i));
    dt=t(2,1)-t(1,1); %it is the same for each case
    a(:,i)=gradient(u(:,i),dt);
end

display(Beta);
display(R);

figure;
subplot(2,1,1);
yyaxis left
plot(t(:,2),a(:,2));
ylabel('acceleration [m/s^2]','FontWeight','bold');
hold on;
yyaxis right
plot(t(:,2),u(:,2));
ylim([-1.1 1.1]);
ylabel('velocity [m/s]','FontWeight','bold');
xlabel('time [s]','FontWeight','bold');
title('Velocity and acceleration for case 2','FontWeight','bold');

subplot(2,1,2);
yyaxis left
plot(t(:,3),a(:,3));
ylabel('acceleration [m/s^2]','FontWeight','bold');
ylim([-1.6 1.6]);
hold on;
yyaxis right
plot(t(:,3),u(:,3));
title('Velocity and acceleration for case 3','FontWeight','bold');
xlabel('time [s]','FontWeight','bold');
ylabel('velocity [m/s]','FontWeight','bold');

%% 9.2 Application to the Egmond fieldwork

%% Sediment transport by waves only

% Sediment transport at each location along the Egmond profile using the SANTOSS model

prof=load('prof1018.txt');
MeanWaterDepth=load('MeanWaterDepth.txt');
midtide=load('midTide.txt');

% Definition of cross-shore coordinates (m)
x = prof(:,1);  

% Definition of zb, bed elevation relative to mean water level (m)
zb = prof(:,2);   

% Definition of the array profile, input argument for BJmodel
profile = [x zb];

% Offshore wave conditions for low, mid, high tide respectively
wlong=-3.27; % alongshore component of the wind velocity, const along the cross-shore direction
wcross=8.55; % cross-shore component of the wind velocity, const along the cross-shore direction
ka=0.022; % apparent bed roughness
nu=0.5; % large-scale mixing coefficient
dzetady=-1.76e-5; % tidally induced 10-100 km scale alongshore slope of the mean sea-surface, const along the cross-shore direction
H130 =2.25; % Significant wave height (m), needed to calculate Hrms0
Hrms0=H130/sqrt(2); % computing Hrms_0 from H13_0
theta0 =39; % Angle of incidence (degrees)
T0 =6.69; % Characteristic period (s)
Zeta =0.09; % Mean water level (m)

% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)

%------------------------------------------------
%       Computation of the alongshore current
%------------------------------------------------

waves = BJmodel(Hrms0,T0,Zeta,theta0,profile,hmin); % computing variables needed for further computation using Battjes and Janssesn model
%% bb
%%% Input parameters for the SANTOSS model

% Wave characteristics
Uw = 1.2;   % orbital velocity amplitude   [m/s]
T = 7;      % period [s]

% Wave shape parameters 
r = 0:0.15:0.6;  % define "r" vector from 0 to 0.6 with step size of 0.15
Nr = length(r);  % number of elements of r
PHI = -pi/2;     % skewed wave 

% Sediment characteristics
D50 = 0.225;  % D50 in mm
D90 = D50;  % D90 in mm
Rhos = 2650; % sediment density in kg/m^3

%%% Initialisation of the output vectors
Qsx = zeros(1,Nr);  % net sediment in the x-direction (m2/s)
Qsy = zeros(1,Nr);  % net sediment in the y-direction (m2/s)
Occ = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the crest
Oct = zeros(1,Nr);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otc = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ott = zeros(1,Nr);  % dimensionless sed. load entrained during the trough period and transported during the trough

for rI = 1:Nr     % loop on the different values of r considered
        
        %1- computation of the time-series of orbital velocity 
        [u(:,rI),t(:,1)]=waveshape(r(rI),PHI,Uw,T);
        
        % 2- computation of the velocity skewness R and the acceleration skewness beta
        [R(rI),Beta(rI)]=velocity_skewness_asymmetry(u(:,rI),t(:,1));

        
        % 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
        Urms(rI)=std(u(:,rI))*100; % multiplied by 100 to get cm/s
        
        % 4- sediment transport calculation
        [Qsx(rI) Qsy(rI) Occ(rI) Oct(rI) Ott(rI) Otc(rI)] = SANTOSSmodel(D50,D90,Rhos,T,Urms(rI),R(rI),Beta(rI),0,0);
end;
    
   