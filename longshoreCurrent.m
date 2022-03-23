function v = longshoreCurrent(profile,dzetady,wcross,wlong,c,theta,Dr,h_tot,st,ka,nu)

% function v = longshoreCurrent(profile,dzetady,wcross,wlong,c,theta,Dr,h_tot,st,ka,nu);
%
% This function computes the cross-shore distribution of the depth-averaged longshore current. Forcing
% include waves, wind and alongshore pressure gradients (due to tides). The function is based on the depth-integrated 
% and longshore-averaged alongshore momentum equation, assuming steady state, no alongshore variations (except that
% of the large scale water level), ..., so that the equation results in an often-used simple 1D equation, in
% which the forcing (wind, wave, waterlevel gradient) balances friction and mixing.
%
% INPUT
%   profile is the used profile (x,z)
%   dzetady is the longshore gradient in waterlevel
%   wcross is the cross-shore wind velocity
%   wlong is the longshore wind velocity
%   c is the wave celerity across the profile
%   theta is the angle of incidence (in degrees) over the profile
%   Dr is the roller dissipation
%   h_tot is the total waterdepth
%	st is the total orbital velocity standard deviation
%   ka is the apparent bed roughness
%   nu is a large-scale mixing coefficient
% OUTPUT
%   v is the cross-shore distribution of the longshore current
%
% This scheme is based on current.m by Falk Feddersen (SIO). It uses the WT83 formulation with alpha = 1.16 to solve
% for the alongshore current. The onshore/offshore boundary condition is dv/dx = 0.
% Gerben Ruessink, May 2000


% --------------------------------------------------------------------
%                           Initialisation
% --------------------------------------------------------------------
x = profile(:,1); % cross-shore coordinates (m)
dx = x(2) - x(1); % grid step (m)
Nx = length(x);   % Amount of cross-shore points to be considered

% endpoint of flow domain
j = find(isnan(Dr));
N = j(1)-1;  % N is the last non-NaN element of the vector Dr

% constants
rho = 1025;  % water density kg/m3
g = 9.81;    % acceleration of gravity m/s2
rhoa = 1.25; % Air density kg/m3 (for computation of the wind stress)
Cd = 0.002;  % Friction factor for the wind

% other usefull stuff
theta = theta*pi/180;                                % Angle of incidence in rad
wtot = sqrt(wcross^2+wlong^2);                       % Total wind velocity 
fc = 0.015*(ka./h_tot).^(1/3);                       % Manning-Strickler formulation


% --------------------------------------------------------------------
%       Computation of the cross-shore vectors of forcing terms
% --------------------------------------------------------------------
% 1. tidal forcing (alongshore gradient dzetady)
Ftide = -g*h_tot*dzetady;

% 2. wind (tsy)                                    
tsy=rhoa*Cd*wtot*wlong;
Fwind = tsy/rho; 

% 3. waves
dSdx=-(sin(theta))./c.*Dr;
Fwaves = -(1/rho).*dSdx;

% compute total force F in m^2/s^2
F = Ftide + Fwind + Fwaves;                         % F: right hand side of energy balance

% --------------------------------------------------------------------
%               Calculation of the alongshore velocity
% --------------------------------------------------------------------

% waterdepth times eddy viscosity
vh = h_tot*nu;

% adjust various matrices to flow domain size
fc = fc(2:N-1);
vh = vh(1:N);
F = F(2:N-1);
st = st(2:N-1);

% set-up tridiagonal matrix Gd, used to compute the diffusion term
% defined such as Gd*v is d/dx(nu*h*dv/dx)
avh = 0.5*(vh(1:N-1)+vh(2:N))/(dx^2);
Gd = diag(-(avh(1:N-2)+avh(2:N-1)),0)+diag(avh(2:N-2),1)+diag(avh(2:N-2),-1);
Gd(1,1) = -avh(2);
Gd(N-2,N-2) = -avh(N-2);

% The whole system is then solved using an iterative process

% First estimate of v in flow domain, Vold 
% To compute Vold, we assume strong-current case, that is to say:
%       - fc <|u|v> = fc |v|v
%       - lateral mixing term neglected (d/dx(nu h dv/dx))
v2 = F./fc;
Vold = sign(v2).*sqrt(abs(v2));

B = wt83(st,Vold);
dB = dwt83(st,Vold);

e = 1;                    % initialise rms error
i = 1;                    % initialise number of iterations
while (i<15) && (e>1e-12) % iterate either 15 times or until rms dynamical error e<1e-12
    a = fc.*(B - dB.*Vold);
    Fd = F - a;
    G  = diag(fc.*dB,0) - Gd;
    Vt = G\Fd;            % compute new estimation of v
    B = wt83(st,Vt);
    dB = dwt83(st,Vt);
    err = F - fc.*B + Gd*Vt;  
    e = sqrt(err'*err/(N-2)); % computate rms error
    Vold = Vt;             % update velocity
    i = i+1;               % update counter iterations
end

% create a vector of size Nx containing the alongshore velocity v
v = NaN(Nx,1);
v(1:N) = [Vt(1);Vt(1:end);Vt(end)];

