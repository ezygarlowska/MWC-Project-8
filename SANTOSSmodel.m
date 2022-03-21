function [Qsx, Qsy, Occ, Oct, Ott, Otc, Rh] = SANTOSSmodel(D50,D90,Rhos,t,Urms,R,beta,ripples,Ux)
% computation of the net sediment transport based on the Santoss practical
% sand transport model
% inputs 
%   D50 [mm]
%   D90 [mm]
%   Rhos sand density [kg/m3]
%   t  wave period [s]
%   Urms root mean square orbital velocity [cm/s]
%   R    degree of velocity skewness [-]
%   beta degree of acceleration skewness [-]
%   ripples   
%         if ripples = 0 no ripples
%         if ripples = 1 ripple characteristics are computed using the method of O'Donoghue et al., 2006, Coastal Engineering.
%   Ux   magnitude of the undertow [cm/s]  (Ux should be >0)
%
% outputs
%   Qsx volumetric net transport  per unit width in x-direction [m2/s]
%   Qsy volumetric net transport  per unit width in y-direction [m2/s]
%   Occ sand load entrained during wave crest period and transported during
%       crest period [-]
%   Oct sand load entrained during wave crest period and transported during
%       trough period [-]
%   Ott sand load entrained during wave trough period and transported during
%       trough period [-]
%   Otc sand load entrained during wave trough period and transported during
%       crest period [-]
%   Rh  ripple height [m] (0 if ripples = 0)



% For this practical, we only consider undertow-type currents 
% and therefore Ang = 180degrees, i.e. opposed to wave propagation
N = length(D50);    
Unet = Ux*ones(N,1);                 % magnitude net current velocity at reference level zref [cm/s]
Ang = 180*ones(N,1);                 % angle between net cur. dir. and positive (onshore) flow dir. [deg]
Zref     = 20*ones(N,1);             % reference level net current velocity [cm]
Tflow(1:N,1)  = {'wc'}  ;            % type of flow 
                                     % w:wave only (=wave tunnel condition), wc: waves+current, sw:surface waves, swc: surface waves+current
% Parameters not-used in our conditions because surface effects are ignored (Tflow=wc) -> defined as NaNs 
hw  = NaN*ones(N,1);                 % wave height [cm]
H =  NaN*ones(N,1);                  % water depth [cm]  

% Ripples:
Rh = zeros(N,1)    ;                 % ripple height [m]
Rl = zeros(N,1)    ;                 % ripple length [m]
if ripples==0      % no ripples
    rip = {'meas'}; % we use the ripples characteristics defined earlier (Rh=0 and Rl=0)
elseif ripples==1  % ripples accounted for
    rip={'comp'};  % computed using adjusted method of O'Donoghue et al., 2006, Coastal Engineering.
end

%--------------------------------------------------------------------------
%%% THE SANTOSS PRACTICAL SAND TRANSPORT MODEL, VERSION 2.07
%--------------------------------------------------------------------------
% will be written down in version 2 paper
% 1 December 2008
% Jebbe van der Werf
% 17 Dec.  2008 version 1C01 (see logboek.doc 17 dec. 2008)
% René Buijsrogge - Add 'calc' (calculation only) option. Read data
%                                         from file. Version santoss1c01
% 19 Jan. 2009 version 1C02
% René Buijsrogge  -  Unet part mayby skipped (see logboek.doc 19 jan 2009)
% 23 Feb. 2009 version 1C03
% René Buijsrogge  - Alternative way Tc and Tt; use sine function
% 25 Feb. 2009 version 1C04
% René Buijsrogge  - Phase I: Load*snelheid concept (alleen aanpassing
% van santoss1c_core.m
% 5 Maart 2009 version 1C04
% René Buijsrogge  - Add Parameter Analysis function (parameter_ana.m)
% 15 June 2009, version 2.00
% René Buijsrogge - new base version
% 18 aug 2009, version 2.02 René Buijsrogge: implement 'current alone'
% 22 sept 2009, version 2.03 René Buijsrogge: implement modifications Dominic van der A:
% Ripple length in ripple predictor (in rip='comp')+ change alpha_r value;
% 15 oct 2009, version 2.04 René Buijsrogge
% Update choice current alone: option No distinction (C3=0) and All-current alone (C3=6)
% Modify transport formula: Add, H and hw to santoss_core + convert hw to m
% 21 october 2009, v2.05.000 René Buijsrogge - Finetuning GWK Schretlen
% 25 nov 2009, version 2.06 René Buijsrogge  - combined mobile and ripple
% roughness
% - ripple predictor modified (smoothened)
% 16 dec 2009, version 2.06.001 Jan Ribberink - representative rms velocity
% instead of peak velocity input for transport model
% 9 jan 2014 Simplifications for wave course, M.T.
%--------------------------------------------------------------------------
                                     
% % For simplicity, we consider that:
warning off all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTION OF SUBMODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bss={'R98'};%bed shear stress model
%R98    Ribberink, 1998, Coastal Engineering.
plag={'R08'};%phase lag model for rippled bed cases
%R08    Ribberink et al., 2008, Journal of Turbulence.
sflt={'D99'};%model for sheet-flow layer thickness
%D99    adjusted method of Dohmen-Janssen, 1999, Ph.D. thesis.
%R08    Ribberink et al., 2008, Journal of Turbulence.
w_c={'W'};%sheet-flow layer based on wave-alone (W) or wave-current stress (WC)
sett=[bss;sflt;w_c;rip;plag];%combined settings


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pcr=1;      % critical phase lag
g=9.81;     % accelaration due to gravity [m/s2]
vis=10^-6;  % kinematic viscosity [m2/s]
delta=0.2;  % reference height for velocities [m]
rhow=1000;  % water density [kg/m3]

% Calibrated coefficients for the model
n=1.20;     % power to the Shields number
m=9.481;    % multiplication factor
alphas=8.01;% phase lag coefficient sheet-flow cases
alphar=9.28;% phase lag coefficient ripple cases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOME PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D50=10^-3*D50;      %from mm --> m
D90=10^-3*D90;      %from mm --> m
H=10^-2*H;          %from cm --> m
Unet=10^-2*Unet;    %from cm/s --> m/s
Urms=10^-2*Urms;    %from cm/s --> m/s
Zref=10^-2*Zref;    %from cm --> m
hw=10^-2*hw;        %from cm --> m
%put NaN Ang's to zero
Ang(isnan(Ang))=0; %[degrees]

%crest and trough orbital velocity assuming second-order Stokes [m/s]
for i=1:length(R)
        u1(i,:)=(2*Urms(i)^2/((2*R(i)-1)^2+1))^0.5;
        u2(i,:)=(2*R(i)-1)*u1(i);
        uwc(i,:)=u1(i)+u2(i);
        uwt(i,:)=u1(i)-u2(i);
        uwcrepr(i,:)=uwc(i,:)*0.5*sqrt(2);
        uwtrepr(i,:)=uwt(i,:)*0.5*sqrt(2);
end

% %wave-induced current not considered as real current [m/s]
%  NaN --> zero value
for i=1:length(Tflow)
    if strcmp(Tflow(i),'w')==1 ||   strcmp(Tflow(i),'sw')==1 ||  isnan(Unet(i))==1;
        Unet(i,:)=0;
    elseif strcmp(Tflow(i),'wc')==1  ||  strcmp(Tflow(i),'swc')==1 || strcmp(Tflow(i),'c')==1;
        Unet(i,:)=Unet(i);
    end
end
%representative orbital velocity [m/s]
for i=1:length(uwc)
        uw(i,:)=(0.5*uwc(i)^2+0.5*uwt(i)^2)^0.5;
        %equivalent to 2^0.5*Urms assuming second-order Stokes with umax=uwc and umin=uwt
end
aw=uw.*t/(2*pi);        % representative orbital velocity amplitude [m]
Delta=(Rhos-rhow)/rhow; % relative sand density [-]
%NaN ripple dimensions --> zero value
Rh(isnan(Rh))=0;  % ripple height
Rl(isnan(Rl))=0;  % ripple length


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE HEIGHT [m], RIPPLE LENGTH [m] (O'Donoghue ea, 2006, C.Eng.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(sett(4),'comp')==1
    %mobility number based on maximum W velocity [-]
    psimax=max(uwc,uwt).^2./(Delta*g.*D50);
    for i=1:length(psimax)
        if D50(i)<=0.00022
            mn(i,:)=0.55;
            ml(i,:)=0.73;
        elseif D50(i)>=0.0003
            mn(i,:)=1.0;
            ml(i,:)=1.0;
        else
            mn(i,:)=0.55+0.45*(10^3*D50(i)-0.22)/(0.3-0.22);
            ml(i,:)=0.73+0.27*(10^3*D50(i)-0.22)/(0.3-0.22);
        end
        if psimax(i)<=190
            nn(i,:)=1.0;
            nl(i,:)=1.0;
        elseif psimax(i)>=240
            nn(i,:)=0.0;
            nl(i,:)=0.0;
        else
            nn(i,:) = 0.5*(1.+cos(2*pi*(psimax(i)-190)/(2*(240-190))));
            nl(i,:) = 0.5*(1.+cos(2*pi*(psimax(i)-190)/(2*(240-190))));
        end
    end
    Rh=aw.*mn.*nn.*(0.275-0.022*psimax.^0.42);
    Rl=aw.*ml.*nl.*(1.970-0.440*psimax.^0.21);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FRICTION FACTOR [-], SHIELDS NUMBERS [-] AND CURRENT VELOCITY AT REFERENCE LEVEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fc,fw,fcw,unet,uc,ut,Sc,St,Swc,Swt,Tc,Tt,Tcu,Tcd,Ttu,Ttd,Scx,Scy,Stx,Sty]=bssR98(N,Ang,aw,D50,D90,delta,Delta,g,H,Unet,uw,uwc,uwt,Zref,Tflow,beta,R,t,rhow,Rhos,Rh,Rl);

[fc,fw,fcw,unet,ucrepr,utrepr,Screpr,Strepr,Swcrepr,Swtrepr,dum1,dum2,dum3,dum4,dum5,dum6,Scxrepr,Scyrepr,Stxrepr,Styrepr]=bssR98(N,Ang,aw,D50,D90,delta,Delta,g,H,Unet,uw,uwcrepr,uwtrepr,Zref,Tflow,beta,R,t,rhow,Rhos,Rh,Rl);

ucxrepr=uwcrepr+unet.*cos(Ang/180*pi);          %x-component ucrest [m/s]
ucyrepr=unet.*sin(Ang/180*pi);                  %y-component ucrest [m/s]
utxrepr=-uwtrepr+unet.*cos(Ang/180*pi);         %x-component utrough [m/s]
utyrepr=unet.*sin(Ang/180*pi);                  %x-component utrough [m/s]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CRITICAL BED SHEAR STRESS [-] (Soulsby, 1997, Dynamics of Marine Sands)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dstar=(Delta*g./vis^2).^(1/3).*D50;             %dimensionless grain-size [-]
Scr=0.3./(1+1.2*Dstar)+0.055*(1-exp(-0.02*Dstar));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SHEET-FLOW LAYER THICKNESS [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(sflt,'D99')==1
    [SFLTc,SFLTt]=sfltD99(D50,Sc,St,Swc,Swt,sett);
elseif strcmp(sflt,'R08')==1
    if strcmp(sett(3),'WC')==1
        SFLTc=D50*10.6.*Sc;
        SFLTt=D50*10.6.*St;
    elseif strcmp(sett(3),'W')==1
        SFLTc=D50*10.6.*Swc;
        SFLTt=D50*10.6.*Swt;
    end
end
SFLT=max(SFLTc,SFLTt);%maximum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FALL VELOCITY SUSPENDED SAND [m/s] (Soulsby, 1997, Dynamics of Marine Sands)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D50s=0.8*D50;%representative grain size [m]
DstarS=(Delta*g./vis^2).^(1/3).*D50s;%dimensionless grain size [-]
wss=vis./D50s.*((10.36^2+1.049*DstarS.^3).^0.5-10.36);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTATION OF SEDIMENT TRANSPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inn_rs=[SFLTc,wss,Tc,SFLTt,Tt,Rh,Screpr,Scr,Strepr,Delta,D50,t,zeros(N,1),ucxrepr,ucyrepr,ucrepr,utxrepr,utyrepr,utrepr,Tcu,Tcd,Ttu,Ttd,Scxrepr,Scyrepr,Stxrepr,Styrepr];%both ripple and sheet-flow cases

[Qsx,Qsy,Occ,Oct,Ott,Otc]=santoss_core(inn_rs,n,alphas,alphar,m,Pcr,Tflow,H,hw,R);

