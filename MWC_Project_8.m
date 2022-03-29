clear; clc; close all;

%% Chapter 9 - Sediment transport modelling

%% 9.1 Preliminary computations- orbital velocity time-series

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
savefig('Matlab8_1');




%% 9.2 Application to the Egmond fieldwork

%% 9.2.1 Sediment transport by waves only

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

% Loading data from the BJ model to calculate 'u'
waves=load('waves.mat');
k = waves.waves(2).k;
eta=waves.waves(2).eta; 
ht=waves.waves(2).ht; 
Hrms=waves.waves(2).Hrms; 
N_last = find(~isnan(eta),1,'last'); 
k=k(1:N_last);
eta=eta(1:N_last);
ht=ht(1:N_last);
Hrms=Hrms(1:N_last);

u=zeros(1000,length(k)); 
for i=(1:length(k))
% Ursell number
Ur_BJ = Ursell(k(i),ht(i),Hrms(i));
% Uw
Uw = Uw_fun(ht(i),Hrms(i),T0);
% Empirical Sk and As
[Sk_BJ_E,As_BJ_E] = empirical_fun(Ur_BJ); 

% We can now compute the orbital velocity 'u'. 
r = computation_r(Sk_BJ_E,As_BJ_E);
phi = computation_phi(Sk_BJ_E,As_BJ_E);
[u(:,i),t] = waveshape(r,phi,Uw,T0); 
%u(:,i) = U; %timeseries for each position 
end 
%% Compute sediment transport at each location 
%%% Input parameters for the SANTOSS model

% Wave characteristics
%Uw = 1.2;   % orbital velocity amplitude   [m/s]
%T = 7;      % period [s]

% Wave shape parameters 
Nu = length(u(1,:));  % number of elements of u
%PHI = phi;     % skewed wave 

% Sediment characteristics
D50 = 0.225;  % D50 in mm
D90 = D50;  % D90 in mm
Rhos = 2650; % sediment density in kg/m^3

%%% Initialisation of the output vectors
Qsx = zeros(1,Nu);  % net sediment in the x-direction (m2/s)
Qsy = zeros(1,Nu);  % net sediment in the y-direction (m2/s)
Occ = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the crest
Oct = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otc = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ott = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the trough

for i=(1:Nu)
[R,Beta]=velocity_skewness_asymmetry(u(:,i),t);
% 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
Urms=std(u(:,i))*100; % multiplied by 100 to get cm/s
% 4- sediment transport calculation
[Qsx(i) Qsy(i) Occ(i) Oct(i) Ott(i) Otc(i)] = SANTOSSmodel(D50,D90,Rhos,T0,Urms,R,Beta,0,0); 
end 

%Evolution of Qsx as a function of x 
figure;
plot(x(1:2472),Qsx); 
ylabel('Q_{sx} [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
xlim([3000 5000]);
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
savefig('Matlab8_2');

%% Sediment transport gradient 
epsilon = 0.6
[QX] = gradient(Qsx);

%Evolution of Qsx as a function of x 
figure;
subplot(3,1,1);
plot(x(1:2472),Qsx); 
xlim([3500 5000]);
ylabel('Q_{sx} [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');

subplot(3,1,2);
plot(x(1:2472),QX); 
xlim([3500 5000]);
ylabel('dQ_{sx}/{dx} [m/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Positional gradient of the volumetric net transport per unit width in x-direction','FontWeight','bold');

subplot(3,1,3);
plot(x(1:2472),zb(1:2472)); 
xlim([3500 5000]);
ylabel('z [m]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Seabed evolution in the x-direction','FontWeight','bold');
savefig('Matlab8_3');

%% Compute the bed level change after 10 days
    %We assume that QX stays constant over time
deltat=3600*24*10; %10 days  
deltazb=-epsilon*QX*deltat;
deltazb=deltazb';
zb10=zb(1:2472)+deltazb(1:2472);

figure; 
plot(x(1:2472),zb(1:2472));
hold on;
plot(x(1:2472),zb10);
ylabel('z [m]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
legend('Z(t=0)','Z(t=10days)');
title('Time evolution of the seabed','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_4');

%% 9.2.2 Sediment transport by waves and current

%Computing undertow 
rho=1000; %Water density [kg/m^3]
E=waves.waves(2).E; %Wave energy
Er=waves.waves(2).Er; %Wave energy roller
c=waves.waves(2).c; %Wave celerity 
Ux=undertow(E(1:2472),Er(1:2472),c(1:2472),ht(1:2472),rho,Hrms(1:2472));
Ux=Ux*100; %in cm/s

%Sediment transport with nonlinearity
Qsxn = zeros(1,Nu);  % net sediment in the x-direction (m2/s)
Qsyn = zeros(1,Nu);  % net sediment in the y-direction (m2/s)
Occn = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the crest
Octn = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otcn = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ottn = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the trough

for i=(1:Nu)
[R,Beta]=velocity_skewness_asymmetry(u(:,i),t);
% 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
Urms=std(u(:,i))*100; % multiplied by 100 to get cm/s
% 4- sediment transport calculation
[Qsxn(i) Qsyn(i) Occn(i) Octn(i) Ottn(i) Otcn(i)] = SANTOSSmodel(D50,D90,Rhos,T0,Urms,R,Beta,0,Ux); 
end 

figure; 
subplot(4,1,1);
plot(x(1:2472),Qsx);
hold on;
plot(x(1:2472),Qsxn);
ylabel('Qsx [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Ux=0','Ux');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,2);
plot(x(1:2472),Qsy);
hold on;
plot(x(1:2472),Qsyn);
ylabel('Qsy [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Ux=0','Ux');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,3);
plot(x(1:2472),Occ);
hold on;
plot(x(1:2472),Occn);
ylabel('Occ','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Ux=0','Ux');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,4);
plot(x(1:2472),Ott);
hold on;
plot(x(1:2472),Ottn);
ylabel('Ott','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Ux=0','Ux');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_5');

%% Sediment transport gradient 
epsilon = 0.6
[QXn] = gradient(Qsxn);

figure; 
plot(x(1:2472),QX);
hold on;
plot(x(1:2472),QXn);
ylabel('X graient of Qsx [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Ux=0','Ux');
title('X graient of Qsx','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_6');

%% Compute the bed level change after 10 days
    %We assume that QX stays constant over time 
deltazbn=-epsilon*QXn*deltat;
deltazbn=deltazbn';
zb10n=zb(1:2472)+deltazbn(1:2472);

figure; 
plot(x(1:2472),zb(1:2472));
hold on;
plot(x(1:2472),zb10);
plot(x(1:2472),zb10n);
ylabel('z [m]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
legend('Z(t=0)','Z(t=10days)','Z(t=10days) with undertow');
title('Time evolution of the seabed','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_7');

%% 9.2.3 Influence of the offshore wave conditions

%% Small waves 
% Loading data from the BJ model to calculate 'u'
Hrms=0.5*Hrms;

u=zeros(1000,length(k)); 
for i=(1:length(k))
% Ursell number
Ur_BJ = Ursell(k(i),ht(i),Hrms(i));
% Uw
Uw = Uw_fun(ht(i),Hrms(i),T0);
% Empirical Sk and As
[Sk_BJ_E,As_BJ_E] = empirical_fun(Ur_BJ); 

% We can now compute the orbital velocity 'u'. 
r = computation_r(Sk_BJ_E,As_BJ_E);
phi = computation_phi(Sk_BJ_E,As_BJ_E);
[u(:,i),t] = waveshape(r,phi,Uw,T0); 
%u(:,i) = U; %timeseries for each position 
end 
%% Sediment transport by waves and current

%Computing undertow 
rho=1000; %Water density [kg/m^3]
E=waves.waves(2).E; %Wave energy
Er=waves.waves(2).Er; %Wave energy roller
c=waves.waves(2).c; %Wave celerity 
Ux=undertow(E(1:2472),Er(1:2472),c(1:2472),ht(1:2472),rho,Hrms(1:2472));
Ux=Ux*100; %in cm/s

%Sediment transport with nonlinearity
Qsxs = zeros(1,Nu);  % net sediment in the x-direction (m2/s)
Qsys = zeros(1,Nu);  % net sediment in the y-direction (m2/s)
Occs = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the crest
Octs = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otcs = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the crest
Otts = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the trough

for i=(1:Nu)
[R,Beta]=velocity_skewness_asymmetry(u(:,i),t);
% 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
Urms=std(u(:,i))*100; % multiplied by 100 to get cm/s
% 4- sediment transport calculation
[Qsxs(i) Qsys(i) Occs(i) Octs(i) Otts(i) Otcs(i)] = SANTOSSmodel(D50,D90,Rhos,T0,Urms,R,Beta,0,Ux); 
end 

figure; 
subplot(4,1,1);
plot(x(1:2472),Qsxs);
hold on;
plot(x(1:2472),Qsxn);
ylabel('Qsx [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','0.5*Hrms');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,2);
plot(x(1:2472),Qsys);
hold on;
plot(x(1:2472),Qsyn);
ylabel('Qsy [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','0.5*Hrms');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,3);
plot(x(1:2472),Occs);
hold on;
plot(x(1:2472),Occn);
ylabel('Occ','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','0.5*Hrms');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,4);
plot(x(1:2472),Otts);
hold on;
plot(x(1:2472),Ottn);
ylabel('Ott','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','0.5*Hrms');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_8');

%% Sediment transport gradient 
epsilon = 0.6
[QXs] = gradient(Qsxs);

figure; 
plot(x(1:2472),QXs);
hold on;
plot(x(1:2472),QXn);
ylabel('X graient of Qsx [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','0.5*Hrms');
title('X graient of Qsx','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_9');

%% Compute the bed level change after 10 days
    %We assume that QX stays constant over time 
deltazbs=-epsilon*QXs*deltat;
deltazbs=deltazbs';
zb10s=zb(1:2472)+deltazbs(1:2472);

figure; 
plot(x(1:2472),zb(1:2472));
hold on;
plot(x(1:2472),zb10s);
plot(x(1:2472),zb10n);
ylabel('z [m]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
legend('Z(t=0)','Z(t=10days)','Z(t=10days) 0.5*Hrms');
title('Time evolution of the seabed','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_10');

%% Large waves 
% Loading data from the BJ model to calculate 'u'
Hrms=4*Hrms; %Twice the original wave height

u=zeros(1000,length(k)); 
for i=(1:length(k))
% Ursell number
Ur_BJ = Ursell(k(i),ht(i),Hrms(i));
% Uw
Uw = Uw_fun(ht(i),Hrms(i),T0);
% Empirical Sk and As
[Sk_BJ_E,As_BJ_E] = empirical_fun(Ur_BJ); 

% We can now compute the orbital velocity 'u'. 
r = computation_r(Sk_BJ_E,As_BJ_E);
phi = computation_phi(Sk_BJ_E,As_BJ_E);
[u(:,i),t] = waveshape(r,phi,Uw,T0); 
%u(:,i) = U; %timeseries for each position 
end 
%% Sediment transport by waves and current

%Computing undertow 
rho=1000; %Water density [kg/m^3]
E=waves.waves(2).E; %Wave energy
Er=waves.waves(2).Er; %Wave energy roller
c=waves.waves(2).c; %Wave celerity 
Ux=undertow(E(1:2472),Er(1:2472),c(1:2472),ht(1:2472),rho,Hrms(1:2472));
Ux=Ux*100; %in cm/s

%Sediment transport with nonlinearity
Qsxl = zeros(1,Nu);  % net sediment in the x-direction (m2/s)
Qsyl = zeros(1,Nu);  % net sediment in the y-direction (m2/s)
Occl = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the crest
Octl = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otcl = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ottl = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the trough

for i=(1:Nu)
[R,Beta]=velocity_skewness_asymmetry(u(:,i),t);
% 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
Urms=std(u(:,i))*100; % multiplied by 100 to get cm/s
% 4- sediment transport calculation
[Qsxl(i) Qsyl(i) Occl(i) Octl(i) Ottl(i) Otcl(i)] = SANTOSSmodel(D50,D90,Rhos,T0,Urms,R,Beta,0,Ux); 
end 

figure; 
subplot(4,1,1);
plot(x(1:2472),Qsxl);
hold on;
plot(x(1:2472),Qsxn);
ylabel('Qsx [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','2*Hrms');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,2);
plot(x(1:2472),Qsyl);
hold on;
plot(x(1:2472),Qsyn);
ylabel('Qsy [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','2*Hrms');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,3);
plot(x(1:2472),Occl);
hold on;
plot(x(1:2472),Occn);
ylabel('Occ','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','2*Hrms');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,4);
plot(x(1:2472),Ottl);
hold on;
plot(x(1:2472),Ottn);
ylabel('Ott','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','2*Hrms');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_11');

%% Sediment transport gradient 
epsilon = 0.6
[QXl] = gradient(Qsxl);

figure; 
plot(x(1:2472),QXl);
hold on;
plot(x(1:2472),QXn);
ylabel('X graient of Qsx [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('Normal conditions','2*Hrms');
title('X graient of Qsx','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_12');

%% Compute the bed level change after 10 days
    %We assume that QX stays constant over time 
deltazbl=-epsilon*QXs*deltat;
deltazbl=deltazbl';
zb10l=zb(1:2472)+deltazbl(1:2472);

figure; 
plot(x(1:2472),zb(1:2472));
hold on;
plot(x(1:2472),zb10l);
plot(x(1:2472),zb10n);
ylabel('z [m]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
legend('Z(t=0)','Z(t=10days)','Z(t=10days) 2*Hrms');
title('Time evolution of the seabed','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_13');

%% 9.2.4 Influence of ripples

%Recalculate u for normal wave conditions
Hrms=0.5*Hrms; %Original wave height

u=zeros(1000,length(k)); 
for i=(1:length(k))
% Ursell number
Ur_BJ = Ursell(k(i),ht(i),Hrms(i));
% Uw
Uw = Uw_fun(ht(i),Hrms(i),T0);
% Empirical Sk and As
[Sk_BJ_E,As_BJ_E] = empirical_fun(Ur_BJ); 

% We can now compute the orbital velocity 'u'. 
r = computation_r(Sk_BJ_E,As_BJ_E);
phi = computation_phi(Sk_BJ_E,As_BJ_E);
[u(:,i),t] = waveshape(r,phi,Uw,T0); 
%u(:,i) = U; %timeseries for each position 
end 
%% Sediment transport by waves and current

%Sediment transport with nonlinearity
Qsxr = zeros(1,Nu);  % net sediment in the x-direction (m2/s)
Qsyr = zeros(1,Nu);  % net sediment in the y-direction (m2/s)
Occr = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the crest
Octr = zeros(1,Nu);  % dimensionless sed. load entrained during the crest period and transported during the  trough
Otcr = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the crest
Ottr = zeros(1,Nu);  % dimensionless sed. load entrained during the trough period and transported during the trough

for i=(1:Nu)
[R,Beta]=velocity_skewness_asymmetry(u(:,i),t);
% 3- computation of the root-mean squared orbital velocity Urms: should be in CM/S!
Urms=std(u(:,i))*100; % multiplied by 100 to get cm/s
% 4- sediment transport calculation
[Qsxr(i) Qsyr(i) Occr(i) Octr(i) Ottr(i) Otcr(i)] = SANTOSSmodel(D50,D90,Rhos,T0,Urms,R,Beta,1,0); 
end 

figure; 
subplot(4,1,1);
plot(x(1:2472),Qsx);
hold on;
plot(x(1:2472),Qsxr);
ylabel('Qsx [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('No ripples','Ripples');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,2);
plot(x(1:2472),Qsy);
hold on;
plot(x(1:2472),Qsyr);
ylabel('Qsy [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('No ripples','Ripples');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,3);
plot(x(1:2472),Occ);
hold on;
plot(x(1:2472),Occr);
ylabel('Occ','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('No ripples','Ripples');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);

subplot(4,1,4);
plot(x(1:2472),Ott);
hold on;
plot(x(1:2472),Ottr);
ylabel('Ott','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('No ripples','Ripples');
title('Effect on undertow on sediment transport','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_14');

%% Sediment transport gradient 
epsilon = 0.6
[QXr] = gradient(Qsxr);

figure; 
plot(x(1:2472),QX);
hold on;
plot(x(1:2472),QXr);
ylabel('X graient of Qsx [m^2/s]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
title('Volumetric net transport per unit width in x-direction','FontWeight','bold');
legend('No ripples','Ripples');
title('X graient of Qsx','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_15');

%% Compute the bed level change after 10 days
    %We assume that QX stays constant over time 
deltazbr=-epsilon*QXn*deltat;
deltazbr=deltazbr';
zb10r=zb(1:2472)+deltazbr(1:2472);

figure; 
plot(x(1:2472),zb(1:2472));
hold on;
plot(x(1:2472),zb10);
plot(x(1:2472),zb10r);
ylabel('z [m]','FontWeight','bold');
xlabel('x [m]','FontWeight','bold');
legend('Z(t=0)','Z(t=10days)','Z(t=10days) with ripples');
title('Time evolution of the seabed','FontWeight','bold');
xlim([4000 5000]);
savefig('Matlab8_16');