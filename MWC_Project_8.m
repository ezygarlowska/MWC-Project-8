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

%% sediment transport analysis
    
   