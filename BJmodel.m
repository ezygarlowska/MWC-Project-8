function waves = BJmodel(Hrms0,T0,Zeta,theta0,profile,hmin)
% BJmodel compute waves over profile
%   WAVES = BJmodel(Hrms0,T0,Zeta,theta0,profile,hmin)
%   returns a structure WAVES with the values of the wave properties at
%   each cross-shore location.
%   BJmodel is based on the Battjes-Jansen model (e.g. Reniers
%   (1997) in Coastal Engineering). 
%
%   Inputs: 
%       Hrms0   root-mean-square wave height at the most offshore location(in m)
%       T0      wave period at the most offshore location (in s)
%       theta0  wave angle to shore normal at the most offshore location (in degrees)
%       Zeta    water level with respect to Mean Sea Level (in m), (+ = above, - = below).
%       profile a [ngridpoints * 2] matrix, with the first column being the
%               cross-shore coordinate (in meters), positive in the onshore
%               direction, and the second column the vertical coordinate (in
%               m), positive above Mean Sea Level.
%
%       hmin    Minimum water depth below which we stop the computation (in m)
%
%   Output: structure WAVES which contains the following field with ngridpoints
%   values in each field (at each cross-shore location in the grid):
%   'Hrms':  root-mean-square wave height
%   'E':     wave energy
%   'eta':   set-up
%   'ht':    total water depth
%   'k':     wave number
%   'c':     wave celerity
%   'n':     proportionality factor (= group velocity divided by
%            celerity)
%   'cg':    group velocity
%   'theta': wave angle
%   'Sxx'    cross-shore component of radiation stress
%   'gamma': parameter in wave-breaking formulation
%   'Hmax':  maximum wave height
%   'Qb':    fraction of breaking waves
%   'Er':    roller energy
%   'Dr':    roller dissipation
%   'st':    standard deviation of near-bed orbital velocity
%   'Dbr':   dissipation due to wave breaking
%   'x':     cross-shore grid
%   'z':     bed level (relative to mean sea level; - is below MSL)
%
%   The values in the WAVES structure are computed for each point in
%   the grid, using the values at the previous (more seaward) point.
%   Because the value of ETA that is computed in each step is also
%   used in the same computation as the difference between the
%   previous and new total water depth HT, it has to be computed
%   iteratively:
%   First, the dissipation due to breaking DBR is computed. DBR is then 
%   subtracted from the wave energy E, and E is also adjusted for the 
%   cross-shore directed component of the group speed of the waves CG.
%   Subsequently, the roller dissipation DR (due to interaction with
%   the wave) is computed and subtracted from the roller energy ER.
%   The amount of energy that goes into breaking DW, is transferred
%   (added) to the roller energy ER. ER is then adjusted for the
%   cross-shore directed component of the wave speed C. Finally, from
%   the radiation stress caused by the gradient in wave motion and the
%   interaction between the wave and the roller, the set-up ETA is
%   computed. The entire process is repeated ITERMAX times, or until
%   the newly computed ETA differs no more from the previously
%   computed ETA than DELETA.

%----------------------------------------------
%               INITIALISATION
%----------------------------------------------

% physical constants
rho = 1025;       % density of water (kg/m3)
g = 9.81;         % graviation of the earth (m/s2)
% model parameters
itermax = 25;     % max number of iterations for re-computing eta
epsEta = 1E-4;    % minimum eta difference (m)
beta = 0.05;      % roller slope (see Ruessink et al. 2001)
alpha = 1;        % alpha parameter in the Dbr formulation

% shorter coding:
ngridpoints = size(profile,1);
x = profile(:,1);
z = profile(:,2);

% initialize variable vectors
nandummy = NaN(ngridpoints,1);
E = nandummy; Hrms= nandummy; Dbr = nandummy; Df = nandummy;
c = nandummy; n = nandummy; cg = nandummy; k = nandummy;
ht = nandummy; Er = nandummy; Dr = nandummy; eta = nandummy;
st = nandummy; Qb = nandummy; theta = nandummy; 
Sxx = nandummy; Hmax = nandummy;

% set values for first gridpoint (offshore boundary)
Hrms(1)  = Hrms0;                     % Root mean square height                 
eta(1)   = 0;                         % Set up
ht(1)    = -z(1) + eta(1) + Zeta;     % Total water depth (includes set-up)
% ----------------- TO BE FILLED IN ---------------
E(1)     = 1/8*rho*g*Hrms(1).^2;                       % Wave energy
k(1)     = k_fun(T0,ht(1));                       % Wave number
c(1)     = phase_fun(T0,ht(1));                       % Phase celerity
n(1)     = n_fun(k(1),ht(1));                       % ratio group velocity/phase celerity (used in the computation of GROUPVELOCITY and RADIATIONSTRESS)
cg(1)    = group_fun(T0,ht(1));                       % Group celerity
theta(1) = deg2rad(theta0);                       % Wave direction in RAD
% ------------------------------------------------
Er(1)    = 0;                         % Roller energy. Needed for computation of Sxx(1), and therefore initialized at 0.
Sxx(1)   = radiationStressXX(n(1),theta(1),E(1),Er(1)); % Radiation stress
st(1)    = stdevOrbital(T0,Hrms(1),k(1),ht(1)); % standart deviation of the orbital velocity (linear theory). Needed to compute the alongshore flow.

% Determine the index of the last wet gridpoint considered
id = find((-z + Zeta) >= hmin);
if isempty(id) % all points are too shallow
    lastWet = 0;
else           % some points are deep enough
    lastWet = id(end) - 1;
end
clear id;

% computation of the gamma parameter
% (assumes that the first grid point is in deep water)
gamma = gammaBS(Hrms(1),k(1));

%-------------------------------------------------------
%           COMPUTATION OF WAVE CHARACTERISTICS 
%-------------------------------------------------------

% compute the values in the waves structure step-wise on the cross-shore grid positions, starting at the seaward end (= x(1))

for gid = 1:lastWet       % loop on the cross-shore positions
    
    % this computation is done iteratively because of set-up
    iter = 1;              % count the number of iterations for set-up computation
    dEta = realmax;        % initialize dEta to a large value to start the iteration process
    eta(gid+1) = eta(gid); % set eta (setup) at gid+1 to the current eta value 
    
    while (iter <= itermax) && (dEta > epsEta)
        
        % maximum wave height
        Hmax(gid) = maxWaveHeight(k(gid),gamma,ht(gid));
        
        % fraction of breaking waves
        if gid == 1
            Qb_guess = 0.005;
        else
            Qb_guess = Qb(gid-1);
        end
        Qb(gid) = fracQbClip(Hrms(gid),Hmax(gid),Qb_guess);
        
        % dissipation due to wave breaking
        Dbr(gid) = dissBreakingBJ(alpha,Qb(gid),T0,Hmax(gid),rho,g);
        
        % computation of some variables at gid+1
        dx = x(gid+1) - x(gid);                                 % cross-shore gridsize
        ht(gid+1) = -z(gid+1) + eta(gid+1) + Zeta;              % total water level
        % ---------------- TO BE FILLED IN ---------------
        k(gid+1) = k_fun(T0,ht(gid+1));                                        % wave number
        c(gid+1) = phase_fun(T0,ht(gid+1));                                        % wave celerity
        n(gid+1) = n_fun(k(gid+1),ht(gid+1));                                        % group velocity/phase celerity
        cg(gid+1) = group_fun(T0,ht(gid+1));                                        % group velocity
        theta(gid+1) = asin(sin(theta(gid)*c(gid+1)/c(gid)));   %We use snell's law                             % wave direction
        %
        % %energy balance
        E(gid+1) = (E(gid)*cg(gid)*cos(theta(gid)))/(cg(gid+1)*cos(theta(gid+1))) - (Dbr(gid))*dx/(cg(gid+1)*cos(theta(gid+1)));                                       % wave energy
        %
        % %compute new Hrms from wave energy
        Hrms(gid+1) = sqrt(E(gid+1)*8/(rho*g));  %We can calculate Hrms from energy, formula (5.1) matlab pdf % root mean square wave height
        % ---------------- END PART TO BE FILLED ---------
        
        % roller dissipation
        Dr(gid) = dissRoller(beta,Er(gid),c(gid),g);
        
        % change in roller flux (it decreases due to roller dissipation caused by interaction with the breaking wave Dr, but
        % increases due to wave breaking Dbr; all energy that is lost to the wave in wave breaking is transferred into the roller).
        dEr = -Dr(gid) + Dbr(gid);
        
        % roller balance
        Er(gid+1) = (2*Er(gid)*c(gid)*cos(theta(gid)) + dEr*dx) /  ...
            (2*c(gid+1)*cos(theta(gid+1))); % cross-shore directed component of the wave celerity
        
        % radiation stress Sxx associated with wave and roller energy
        Sxx(gid+1) = radiationStressXX(n(gid+1),theta(gid+1),E(gid+1),Er(gid+1));
        
        % change in Sxx
        dSxxdx = (Sxx(gid+1) - Sxx(gid))/dx;
        
        % set-up associated with radiation stress (wave and roller)
        oldEstimate = eta(gid+1);
        eta(gid+1) = eta(gid) - dSxxdx*dx/(rho*g*mean(ht(gid:gid+1)));
        
        % standard deviation orbital flow
        st(gid+1) = stdevOrbital(T0,Hrms(gid+1),k(gid+1),ht(gid+1));        
    
        % prepare for next iteration
        dEta = abs(eta(gid+1) - oldEstimate);
        iter = iter + 1;  % count iterations
    end
        
    if E(gid+1) <= 0
        disp('We have run out of wave energy. This should not happen with a sufficiently large hmin and a sufficiently dense grid');
        break
    end
end

%-----------------------------------------
%               OUTPUTS
%-----------------------------------------

% create output structure WAVES; convert theta back to degrees
waves = struct('E',E, 'Hrms',Hrms, 'Dbr',Dbr, 'c',c, ...
    'n',n, 'cg',cg, 'k',k, 'ht',ht, 'Er',Er, 'Dr', Dr, 'eta',eta, ...
    'Qb',Qb, 'theta',theta*180/pi, 'gamma',gamma, 'st', st, ...
    'Sxx',Sxx, 'Hmax',Hmax, 'x',x, 'z',z);