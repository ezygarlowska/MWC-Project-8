%--------------------------------------------------------------------------
%   COMPUTATION OF THE BED SHEAR STRESSES USING THE METHOD OF RIBBERINK
%   (1998);  ACCELERATION-SKEWNESS method by van der A (2009)                                                            
%   Created by Jebbe van der Werf on 05-07-2006
%   Modified by Jebbe van der Werf on 28-06-2007  
%   May 2009:   Phase 2 model
%   Modified by René Buijsrogge, June 2009
%  19 aug 2009, version 2-02; R. Buijsrogge - update z0c in 'current alone'  part
%  22 sept 2009, version 2-03 René Buijsrogge: implement modifications Dominic van der A:
%  Use ripple length and ripple heigth
%  20 oct 2009, version 2.04, Dominic van der A/René Buijsrogge
%  - Modify Initial current roughness, mobile bed rougness for sheet-flow
%    conds, additional roughness if wave-aaveraged total Shields > 1 m.
%  - Modify X,awc,wat expression of accelereation skewness
%  - Add effect of wave reynolds stress (bl streaming) on near bed crest and
%    trough periode Tc,Tt,Tcu,Ttu
%  21 october 2009, version 2.05 René Buijsrogge
%  - crest period extension for horizontal particle displacement.
%  - Tuning factor 0.55 from measurements GWK Schretlen 2010
%  25 nov 2009, version 2.06 René Buijsrogge
%  - combined mobile and ripple roughness
%  8 jan 2010, René Buijsrogge
%  -  correction weighing waves and current
%  Copyright 2009 University of Twente, Water Engineering and Management, Enschede, The Netherlands
%  18 jan 2010. Put factor mu now in a parameter name mu in stead of number
%---------------------------------------------------------------------------------------------------

function [fc,fw,fcw,unet,uc,ut,Sc,St,Swc,Swt,Tc,Tt,Tcu,Tcd,Ttu,Ttd,Scx,Scy,Stx,Sty]=bssR98(N,Ang,aw,D50,D90,delta,Delta,g,H,Unet,uw,uwc,uwt,Zref,Tflow,b,r,t,rhow,Rhos,Rh,Rl)

% Initialize formal parameters
fc = zeros(size(Ang));
fw = zeros(size(Ang));
fcw = zeros(size(Ang));
unet = zeros(size(Ang));
uc = zeros(size(Ang));
ut = zeros(size(Ang));
Sc = zeros(size(Ang));
St = zeros(size(Ang));
Swc = zeros(size(Ang));
Swt = zeros(size(Ang));
Tc = zeros(size(Ang));
Tt = zeros(size(Ang));
Tcu = zeros(size(Ang));
Tcd = zeros(size(Ang));
Ttu = zeros(size(Ang));
Ttd = zeros(size(Ang));
Scx = zeros(size(Ang));
Scy = zeros(size(Ang));
Stx = zeros(size(Ang));
Sty = zeros(size(Ang));

% Initialize local variables
mu = 6.; % in combination with deltas sheetflow layer thickness in sfltD99
ksw1 = zeros(size(Ang));
fw1 = zeros(size(Ang));
fc1 = zeros(size(Ang));
dum = zeros(size(Ang));
ksw = zeros(size(Ang));
ksc = zeros(size(Ang));
ksc1 = zeros(size(Ang));
z0c1 = zeros(size(Ang));
z0c = zeros(size(Ang));
fcc = zeros(size(Ang));
theta = zeros(size(Ang));
alpha = zeros(size(Ang));% alpha correction for wave alone or wave + current
alpha_SW = zeros(size(Ang));% alpha correction reynolds stress for surface waves 
theta1 = zeros(size(Ang));
theta2 = zeros(size(Ang));
ustarc = zeros(size(Ang));
X = zeros(size(Ang));
awc = zeros(size(Ang));
awt = zeros(size(Ang));
fwc = zeros(size(Ang));
fwt = zeros(size(Ang));
fcwc = zeros(size(Ang));
fcwt = zeros(size(Ang));

for i=1:N;  % Main loop for all experiments
  if strcmp(Tflow(i),'c')==0  
    %================================================================    
    % (III) Bed roughness (mobile or ripples) and friction coefficient for 
    %             WAVES ALONE or WAVES + CURRENT
    %================================================================  
    %- friction coefficient for waves and for current
    %- mean magnitude of bed shear stress
%     if Rh(i,:)~=0
%         % rippled bed roughness [m]
%         ksw1(i)=D50(i)+Rh(i)*Rh(i)/Rl(i);
%     else
        % initial roughness sheet flow regime [m]
        ksw1(i)=D50(i);
%     end
    %wave friction formula of Swart (1974) [-]
    % for i=1:length(ksw1)
    if ksw1(i)/aw(i) < 0.63
        fw1(i)=exp(5.213*(ksw1(i)/aw(i))^0.194-5.977);
    else
        fw1(i)=0.3;
    end
   %end
    %initial current roughness [m]
    p_corr = 0.4; % Correction factor p. Used for form roughness ripples
    if Rh(i,:)~=0
       ksc1(i)=3*D90(i)+Rh(i)*Rh(i)/Rl(i)*p_corr;
    else
        ksc1(i)=3*D90(i);
    end
    z0c1(i)=ksc1(i)/30; 
    %current fricion factor assuming logarithmic current profile [-]
    % for i=1:length(ksc1)
    if Unet(i)==0
        fc1(i)=0;
    elseif Zref(i)==1000 %depth-averaged current
        fc1(i)=2*(0.4/(log(H(i)/z0c1(i))-1+z0c1(i)/H(i)))^2;
    else
        fc1(i)=2*(0.4/log(Zref(i)/z0c1(i)))^2;
    end
    % end
    %initial wave-averaged total bed shear stress [-] 
    theta1(i)=0.5*fc1(i)*Unet(i)^2./(Delta(i)*g.*D50(i))+0.25*fw1(i)*uw(i)^2./(Delta(i)*g*D50(i));
    % - Mobile bed roughness
    % - Friction coefficient for waves and for current
    % - mean magnitude of bed shear stress
  
    %start of Shields loop
    % for i=1:length(theta1)
%     if Rh(i)==0 % mobile bed roughness only for sheet-flow conds        
       clear dum
       j=1;    
       dum(j)=theta1(i);    
       %additional roughness if wave-averaged total Shields > 1 [m]
       if D50(i)<=0.00015
           ksw(i)=max(ksw1(i),D50(i)*(mu+6*(dum(j)-1)));       
           ksc(i)=max(3*D90(i),D50(i)*(mu+6*(dum(j)-1)));        
       elseif D50(i)>=0.00020
           ksw(i)=max(ksw1(i),D50(i)*(1+6*(dum(j)-1)));       
           ksc(i)=max(3*D90(i),D50(i)*(1+6*(dum(j)-1)));        
       else
           ksw(i)=max(ksw1(i),D50(i)*(mu+(D50(i)-0.00015)*(1-mu)/(0.00020-0.00015)+6*(dum(j)-1)));       
           ksc(i)=max(3*D90(i),D50(i)*(mu+(D50(i)-0.00015)*(1-mu)/(0.00020-0.00015)+6*(dum(j)-1)));  
       end
       z0c(i)=ksc(i)/30;
       if ksw(i)/aw(i) < 0.63
           fw(i)=exp(5.213*(ksw(i)/aw(i))^0.194-5.977);
       else
           fw(i)=0.3;
       end   
       if Unet(i)==0
           fcc(i)=0;
       elseif Zref(i)==1000
           fcc(i)=2*(0.4/(log(H(i)/z0c(i))-1+z0c(i)/H(i)))^2;
       else
           fcc(i)=2*(0.4/log(Zref(i)/z0c(i)))^2;
       end      
       %wave-averaged total bed shear stress [-] 
       theta2(i)=0.5*fcc(i)*Unet(i)^2/(Delta(i)*g*D50(i))+0.25*fw(i)*uw(i)^2/(Delta(i)*g*D50(i));
       j=j+1;
       dum(j)=theta2(i); 
       %while loop with 0.001 difference as stop criterion
       while abs(dum(j)-dum(j-1))>0.001       
            if D50(i)<=0.00015
               ksw(i)=max(ksw1(i),D50(i)*(mu+6*(dum(j)-1)));       
               ksc(i)=max(3*D90(i),D50(i)*(mu+6*(dum(j)-1)));        
            elseif D50(i)>=0.00020
               ksw(i)=max(ksw1(i),D50(i)*(1+6*(dum(j)-1)));       
               ksc(i)=max(3*D90(i),D50(i)*(1+6*(dum(j)-1)));        
            else
               ksw(i)=max(ksw1(i),D50(i)*(mu+(D50(i)-0.00015)*(1-mu)/(0.00020-0.00015)+6*(dum(j)-1)));       
               ksc(i)=max(3*D90(i),D50(i)*(mu+(D50(i)-0.00015)*(1-mu)/(0.00020-0.00015)+6*(dum(j)-1)));        
            end
            z0c(i)=ksc(i)/30;       
            if ksw(i)/aw(i) < 0.63
                fw(i)=exp(5.213*(ksw(i)/aw(i))^0.194-5.977);
            else
                fw(i)=0.3;
            end
            if Unet(i)==0
                fcc(i)=0;
            elseif Zref(i)==1000
                fcc(i)=2*(0.4/(log(H(i)/z0c(i))-1+z0c(i)/H(i)))^2;
            else
                fcc(i)=2*(0.4/log(Zref(i)/z0c(i)))^2;
            end
            j=j+1;
            dum(j)=0.5*fcc(i)*Unet(i)^2/(Delta(i)*g*D50(i))+0.25*fw(i)*uw(i)^2/(Delta(i)*g*D50(i));
       end
       theta(i)=dum(j);%total Shields [-]
       if Rh(i)~=0
         % rippled bed roughness [m]
         ksw(i)=ksw(i)+Rh(i)*Rh(i)/Rl(i)*p_corr;% form roughness ripples with correction factor p
         if ksw(i)/aw(i) < 0.63
           fw(i)=exp(5.213*(ksw(i)/aw(i))^0.194-5.977);
         else
           fw(i)=0.3;
         end  
       end
%     end
    ustarc(i)=(0.5*fcc(i)).^0.5.*Unet(i);      %friction velocity [m/s]
    unet(i)=ustarc(i)/0.4*log(delta/z0c(i));   %net current strength at z=delta [m/s]
    nu_corr=3.0; % Correctiefactor weging golven en stroom
    alpha(i)=nu_corr*unet(i)/(nu_corr*unet(i)+uw(i));     %relative current strength [-]
    % - unet at z=delta (from Unet at z=zref)
    %for i=1:length(fcc)
    if unet(i)==0
          fc(i)  =0;
          fcw(i)=fw(i);%combined wave-current friction factor [-]
          %magnitude total velocity vector at times of wave crest and wave
          %trough [m/s]
          uc(i) =uwc(i);
          ut(i) =uwt(i);    
    else
          %current friction factor corresponding to unet such that bed shear
          %stress stays the same
          fc(i)=fcc(i)*Unet(i)^2/unet(i)^2;
          %combined wave-current friction coefficient fcw using formula
          % Madsen & Grant (1976)
          fcw(i) =alpha(i).*fc(i)+(1-alpha(i)).*fw(i);
          uc(i)   =((uwc(i)+unet(i)*cos(Ang(i)/180*pi))^2+(unet(i)*sin(Ang(i)/180*pi))^2)^0.5;
          ut(i)   =((unet(i)*cos(Ang(i)/180*pi)-uwt(i))^2+(unet(i)*sin(Ang(i)/180*pi))^2)^0.5;   
    end
    % end
    %================================================================  
    % PARTIAL WAVE CREST AND TROUGH PERIODS [s]
    %================================================================  
    [Tc(i),Tt(i),Tcu(i),Tcd(i),Ttu(i),Ttd(i)]=wctp(Ang(i),r(i),t(i),unet(i),uwc(i),uwt(i),b(i));
    %================================================================
    % ACCELERATION-SKEWNESS
    %================================================================
    %for i=1:length(b)
    if (b(i) ~= 0.5)
        % - wave friction coefficient for crest and trough 
        % - combined wave-current friction coefficent for crest and trough  
        X(i)=(2*Tcu(i)/Tc(i)); %alternative acceleration skewness paramater [-]        
        awc(i)=X(i).^2*aw(i); %equivalent excursion amplitude for crest [m]
        awt(i)=(2-X(i).^2)*aw(i); %equivalent excursion amplitude for trough [m]
        %% crest friction factor
        if ksw(i)/awc(i) < 0.63
           fwc(i)=exp(5.213*(ksw(i)/awc(i))^0.194-5.977);
        else
           fwc(i)=0.3;
        end
        %% trough friction factor
        if ksw(i)/awt(i) < 0.63
            fwt(i)=exp(5.213*(ksw(i)/awt(i))^0.194-5.977);
        else
            fwt(i)=0.3;
        end
        if unet(i)==0
            fc(i)  =0;
            fcwc(i)=fwc(i);%combined wave-current friction factor crest [-]
            fcwt(i)=fwt(i);%combined wave-current friction factor trough [-]
        else         
            %combined wave-current friction coefficient at crest fcwc and trough fcwt
            %using formula Madsen & Grant (1976)
            fcwc(i) =alpha(i).*fc(i)+(1-alpha(i)).*fwc(i);
            fcwt(i) =alpha(i).*fc(i)+(1-alpha(i)).*fwt(i);
        end
        % Bed shear stress for acceleration-skewed waves with/without 
        % a  current Tc, Tt and x,y components Stx, Sty) 
        % bed shear stresses due to wave alone!
        Swc(i) = 0.5*fwc(i)*uwc(i)^2./(Delta(i)*g*D50(i));
        Swt(i) = 0.5*fwt(i)*uwt(i)^2./(Delta(i)*g*D50(i));
        % total bed shear stresses!
        Sc(i)   = 0.5*fcwc(i)*uc(i)^2./(Delta(i)*g*D50(i));
        St(i)   = 0.5*fcwt(i)*ut(i)^2./(Delta(i)*g*D50(i));
        Scx(i) = Sc(i) * (uwc(i) + unet(i)*cos(Ang(i)/180*pi))/uc(i);
        Scy(i) = Sc(i) * unet(i)*sin(Ang(i)/180*pi)/uc(i);
        Stx(i) = St(i)  * (-uwt(i) + unet(i)*cos(Ang(i)/180*pi))/ut(i);
        Sty(i)=  St(i)  * unet(i)*sin(Ang(i)/180*pi)/ut(i); 
    else
        % b = 0.5
        % (IV) Bed shear stress for velocity-skewed or sine waves with/without a current Tc, Tt and x,y components Stx, Sty
        % bed shear stresses due to wave alone!
        Swc(i) = 0.5*fw(i)*uwc(i)^2./(Delta(i)*g*D50(i));
        Swt(i) = 0.5*fw(i)*uwt(i)^2./(Delta(i)*g*D50(i));
        %total bed shear stresses!
        Sc(i)   = 0.5*fcw(i)*uc(i)^2./(Delta(i)*g*D50(i));
        St(i)   = 0.5*fcw(i)*ut(i)^2./(Delta(i)*g*D50(i));
        Scx(i) = Sc(i) * (uwc(i) + unet(i)*cos(Ang(i)/180*pi))/uc(i);
        Scy(i) = Sc(i) * unet(i)*sin(Ang(i)/180*pi)/uc(i);
        Stx(i) = St(i)  * (-uwt(i) + unet(i)*cos(Ang(i)/180*pi))/ut(i);
        Sty(i)=  St(i)  * unet(i)*sin(Ang(i)/180*pi)/ut(i);
    end
    %end
    %================================================================
    % SURFACE WAVES
    %================================================================
    %
    if (strcmp(Tflow(i),'sw') || strcmp(Tflow(i),'swc'))
        %==============================================================
        % (V) Wave Reynolds stress
        %==============================================================
        % Calculation wave propagation speed c
        KSI = 4*pi^2*H(i)/(g*t(i)^2);  % (H=depth, t=wave period)
        if  KSI <= 1.
            eta = sqrt(KSI) * (1. + 0.2 * KSI);
        elseif KSI > 1.
            eta = KSI * (1.+0.2*exp(2. - 2*KSI));
        else
            disp(['Unknown option. eta= ',eta]);
        end
        L= 2* pi* H(i) / eta;
        c = L / t(i);     
        % Wave Reynolds stress and Shields
        %  - alpha = 0.424 (SINE WAVES)
        %  - friction coefficient = fcw (mobile bed)
        alpha_SW(i) = 0.424;     
        tauwRe = rhow * fcw(i) * alpha_SW(i) * uw(i)^3 / 2/ c;
        SwRe    = tauwRe / ((Rhos(i) - rhow) * g * D50(i));
        % ==============================================================
        % (VI) crest period extension for horizontal particle displacement.
        % Tuning factor 0.55 from measurements GWK Schretlen 2010
        % ==============================================================
        Tcu(i) = Tcu(i)/Tc(i)*(Tc(i) + t(i)*0.55*uw(i)/(c*pi-0.55*2*uw(i)) );
        Tc(i)  = Tc(i) + t(i)*0.55*uw(i)/(c*pi-0.55*2*uw(i)) ;     
        Tcd(i) = Tc(i)-Tcu(i);
        Ttu(i) = Ttu(i)/Tt(i)*(Tt(i) - t(i)*0.55*uw(i)/(c*pi+0.55*2*uw(i)) );
        Tt(i)  = Tt(i) - t(i)*0.55*uw(i)/(c*pi+0.55*2*uw(i)) ;    
        Ttd(i) = Tt(i)-Ttu(i);    

        % ==============================================================
        % (VII) Bed shear stress for surface waves with/without a current
        % Tc, Tt and x,y components Stx, Sty
        % ==============================================================
        Scx(i) = Scx(i) + SwRe;
        Stx(i) = Stx(i) + SwRe;
        Sc(i) = sqrt(Scx(i)^2 + Scy(i)^2);
        St(i) = sqrt(Stx(i)^2 + Sty(i)^2);
    end
else
    % ==============================================================
    % CURRENT ALONE PART I AND II
    % ==============================================================
    disp('Current alone found! ')
    % =================================================
    % (I) Mobile bed roughness and friction coefficient for current alone
    % =================================================
    %initial current roughness [m]
    ksc1(i)=3*D90(i);
    z0c1(i)=ksc1(i)/30;  
    %current fricion factor assuming logarithmic current profile [-]
    %for i=1:length(ksc1)
    if Unet(i)==0
       fc1(i)=0;
    elseif Zref(i)==1000%depth-averaged current
       fc1(i)=2*(0.4/(log(H(i)/z0c1(i))-1+z0c1(i)/H(i)))^2;
    else
       fc1(i)=2*(0.4/log(Zref(i)/z0c1(i)))^2;
    end
    %end
    % initial total bed shear stress [-] 
    theta1(i)  =  0.5*fc1(i)*Unet(i)^2./(Delta(i)*g*D50(i));
    % - Mobile bed roughness
    % - Friction coefficient for waves and for current
    % - mean magnitude of bed shear stress
    clear dum
    j=1;    
    dum(j)=theta1(i); 
    %additional roughness if wave-averaged total Shields > 1 [m]      
    if D50(i)<=0.00015
       ksc(i)=max(3*D90(i),D50(i)*(mu+6*(dum(j)-1)));        
    elseif D50(i)>=0.00020    
       ksc(i)=max(3*D90(i),D50(i)*(1+6*(dum(j)-1)));        
    else   
       ksc(i)=max(3*D90(i),D50(i)*(mu+(D50(i)-0.00015)*(1-mu)/(0.00020-0.00015)+6*(dum(j)-1)));        
     end          
     z0c(i)=ksc(i)/30;       
     if Zref(i)==1000
        fcc(i,:)=2*(0.4/(log(H(i)/z0c(i))-1+z0c(i)/H(i)))^2;
     else
        fcc(i,:)=2*(0.4/log(Zref(i)/z0c(i)))^2;
     end
     % total bed shear stress [-] 
     theta2(i) = 0.5*fcc(i) * Unet(i)^2/(Delta(i)*g*D50(i));
     j=j+1;
     dum(j)=theta2(i);
           
     %while loop with 0.001 difference as stop criterion
     while abs(dum(j)-dum(j-1))>0.001       
         if D50(i)<=0.00015
            ksc(i)=max(3*D90(i),D50(i)*(mu+6*(dum(j)-1)));        
         elseif D50(i)>=0.00020    
            ksc(i)=max(3*D90(i),D50(i)*(1+6*(dum(j)-1)));        
         else   
            ksc(i)=max(3*D90(i),D50(i)*(mu+(D50(i)-0.00015)*(1-mu)/(0.00020-0.00015)+6*(dum(j)-1)));        
         end            
         z0c(i)=ksc(i)/30;
         if Zref(i)==1000
            fcc(i)=2*(0.4/(log(H(i)/z0c(i))-1+z0c(i)/H(i)))^2;
         else
            fcc(i)=2*(0.4/log(Zref(i)/z0c(i)))^2;
         end
         j=j+1;
         dum(j)=0.5*fcc(i)*Unet(i)^2/(Delta(i)*g*D50(i));
     end
     theta(i)=dum(j);%total Shields [-]           
    %====================================================================
    % (II) Bed shear stress for only current: Tc, Tt and x,y components Stx, Sty
    %=====================================================================
    ustarc(i)=(0.5*fcc(i))^0.5*Unet(i);      %friction velocity [m/s]
    unet(i)=ustarc(i)/0.4.*log(delta/z0c(i)); %net current strength at z=delta [m/s]
    % for i=1:length(fcc)
    %current friction factor corresponding to unet such that bed shear
    %stress stays the same
    fc(i)=fcc(i)*Unet(i)^2/unet(i)^2;
    % end
    Sc(i) = 0.5 * fc(i) * unet(i)^2 / (Delta(i) * g * D50(i));
    Scx(i) = Sc(i) * cos(Ang(i)/180/pi);
    Scy(i) = Sc(i) * sin(Ang(i)/180/pi);
    Stx(i) = 0;
    Sty(i) = 0;
    Tc(i) = t(i);
    Tt(i) = 0;
    Tcu(i) = 0;
    Tcd(i) = 0;
    Ttu(i) = 0;
    Ttd(i) = 0;
  end
end  % end of main loop for all experiments