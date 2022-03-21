function [Qsx,Qsy,Occ,Oct,Ott,Otc]=santoss_core(inn,n,alphas,alphar,m,Pcr,Tflow,H,hw,r)
%--------------------------------------------------------------------------
%%% THE SANTOSS PRACTICAL SAND TRANSPORT MODEL, VERSION 2.05         
%--------------------------------------------------------------------------
%%% will be written down in version 2 paper
%%% 1 December 2008
%%% Jebbe van der Werf
%%% 19 December 2008
%%% Modified by René Buijsrogge - Add extra mmode ('calc')  calculations only. Read data
%%%                                         from file. Version santoss1c01
%%% 5 February 2009
%%% René Buijsrogge - Bugfix Otc in case Pt = Inifinite
%%% 25 Feb. 2009 version 1C04
%%% René Buijsrogge  - Phase I: Load*snelheid concept 
%%% 9  June 2009, v2.00.000 René Buijsrogge - new base version
%%% 18 Aug 2009, v2.02.000 René Buijsrogge - implement 'current alone'
%%% 15 oct 2009, v2.04.000  René Buijsrogge/Dominic van der A
% - modify test 'cur alone' for use 'cur alone' and 'All-cur alone' together. 
% - modify 1) acc. skewness load weights
% - modify 2) settling vel. correction for surface waves (1st order)  
%%% 21 october 2009, v2.05.000 René Buijsrogge -  Finetuning GWK
% - settling vel. correction for surface waves using maximum vertical
% - orbital velocity for real waves with 2nd order Stokes theory
%--------------------------------------------------------------------------
warning off all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ASSIGNING INPUT
SFLTc=inn(:,1);
wss=inn(:,2); 
Tc=inn(:,3);
SFLTt=inn(:,4);
Tt=inn(:,5);
Rh=inn(:,6);
Sc=inn(:,7);
Scr=inn(:,8);
St=inn(:,9);
Delta=inn(:,10);
D50=inn(:,11);
t=inn(:,12);
ucx=inn(:,14);
ucy=inn(:,15);
uc=inn(:,16);
utx=inn(:,17);
uty=inn(:,18);
ut=inn(:,19);
Tcu=inn(:,20);
Tcd=inn(:,21);
Ttu=inn(:,22);
Ttd=inn(:,23);
Scx=inn(:,24);
Scy=inn(:,25);
Stx=inn(:,26);
Sty=inn(:,27);
g=9.81;

% Initialize
worbc=zeros(length(Rh),1);
worbt=zeros(length(Rh),1);
Psc=zeros(length(Rh),1);
Pst=zeros(length(Rh),1);
Prc=zeros(length(Rh),1);
Prt=zeros(length(Rh),1);
Pc=zeros(length(Rh),1);
Pt=zeros(length(Rh),1);
X= zeros(length(Rh),1);
Phicx= zeros(length(Rh),1);
Phix= zeros(length(Rh),1);
fac2= zeros(length(Rh),1);
fac5= zeros(length(Rh),1);
fac10= zeros(length(Rh),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATION PHASE-LAG PARAMETERS [-]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Rh)
    %calculate maximum vertical orbital velocity for real waves (zero for other)
    %with 2nd order Stokes theory. Tuning factor 1.5 for transport
    %measurements GWK Schretlen 2010
    if (strcmp(Tflow(i),'sw') || strcmp(Tflow(i),'swc'))
        worbc1= pi*hw(i)*SFLTc(i)/t(i)/H(i);
        worbt1= pi*hw(i)*SFLTt(i)/t(i)/H(i);    
        worbc2= worbc1*2*(r(i)+r(i)-1);  
        worbt2= worbt1*2*(r(i)+r(i)-1);       
        eps_corr=3.0; % correction factor level vertical orbital velocity, calibration param
        worbc(i)= (1/8)*worbc1*sqrt(64-(-worbc1+sqrt(worbc1^2+32*worbc2^2))^2/(worbc2^2))+...
                     worbc2*sin(2*acos((1/8)* (-worbc1+sqrt(worbc1^2+32*worbc2^2))/worbc2));
        worbc(i)=worbc(i)*eps_corr;
        worbt(i)= (1/8)*worbt1*sqrt(64-(-worbt1+sqrt(worbt1^2+32*worbt2^2))^2/(worbt2^2))+...
                     worbt2*sin(2*acos((1/8)* (-worbt1+sqrt(worbt1^2+32*worbt2^2))/worbt2)); 
        worbt(i)= worbt(i)*eps_corr;        
    else
        worbc(i)=0;
        worbt(i)=0;
    end

    Psc(i,:)=alphas*SFLTc(i)/((wss(i)+worbc(i))*2*(Tc(i)-Tcu(i)));
    delta_w=max((wss(i)-worbt(i)),0);% test to prevent negative vertical settling velocity
    Pst(i,:)=alphas*SFLTt(i)/(delta_w*2*(Tt(i)-Ttu(i)));
    Prc(i,:)=alphar*Rh(i)/((wss(i)+worbc(i))*2*(Tc(i)-Tcu(i)));
    Prt(i,:)=alphar*Rh(i)/(delta_w*2*(Tt(i)-Ttu(i)));  
    if Rh(i,:)==0
        Pc(i,:)=Psc(i);
        Pt(i,:)=Pst(i);        
    else
        Pc(i,:)=Prc(i);
        Pt(i,:)=Prt(i);         
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD COMPONENTS [-]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Oc=m*max((Sc-Scr),0).^n;%load entrained during crest period
Occ=min(Pcr./Pc,1).*Oc;%load entrained and transported during crest period
Oct=max((Pc-Pcr)./Pc,0).*Oc;%load entrained during crest period and transported during trough period
Ot=m*max((St-Scr),0).^n;%load entrained during trough period
Ott=min(Pcr./Pt,1).*Ot;%load entrained and transported during trough period    
tmp= zeros(length(Tt),1);
for i=1:length(Tt)
  if  abs(Tt(i,:)) < eps % check zero value Ttrough
      tmp(i) = 1. ;   % limit goes to 1. Normal calculation goes to NaN and Otc becomes 0 and not Ot
  else
      tmp(i) = (Pt(i,:)-Pcr)./Pt(i,:) ;% normal calculation
  end
end
Otc=max(tmp,0).*Ot;%entrained during trough period and transported during crest period  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TRANSPORT COMPONENTS [-]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Rh)
    X(i,:)=2*Tcu(i)/Tc(i); % alternative acceleration skewness paramater [-]
end
for i=1:length(Rh)
  if strcmp(Tflow(i),'c')
      Phicx(i)=Scx(i)/sqrt(Sc(i))*Oc(i);                %dimensionless transport in ucx-direction current alone
  else
      Phicx(i)=Scx(i)/sqrt(Sc(i))*(Occ(i)+1/X(i)*Otc(i)); %dimensionless transport in ucx-direction 
  end
end
Phicy=Scy./sqrt(Sc).*(Occ+1./X.*Otc);         %dimensionless transport in ucy-direction
Phitx=Stx./sqrt(St).*(Ott+1./(2-X).*Oct);     %dimensionless transport in utx-direction
Phity=Sty./sqrt(St).*(Ott+1./(2-X).*Oct);     %dimensionless transport in uty-direction
for i=1:length(Rh)
  if strcmp(Tflow(i),'c')
      Phix(i)= Phicx(i);                       %dimensionless transport in x-direction current alone
  else
      Phix(i)=(Tc(i)*Phicx(i)+Tt(i)*Phitx(i))/t(i);   %dimensionless transport in x-direction
  end
end
Phiy=(Tc.*Phicy+Tt.*Phity)./t;         %dimensionless transport in y-direction
Qsx=Phix.*(Delta.*g.*D50.^3).^0.5;     %transport in x-direction [m2/s]
Qsy=Phiy.*(Delta.*g.*D50.^3).^0.5;     %transport in y-direction [m2/s]

