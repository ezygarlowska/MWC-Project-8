%--------------------------------------------------------------------------
%%% Santoss Transport Model
%--------------------------------------------------------------------------
function [Tc,Tt,Tcu,Tcd,Ttu,Ttd]=wctp(Ang,r,t,unet,uwc,uwt,b)
%--------------------------------------------------------------------------
%%% PARTIAL WAVE CREST AND TROUGH PERIODS [s]
%--------------------------------------------------------------------------
%%% 19 Mei 2009, version 2-00
%%% René Buijsrogge - new base version.
%%% Move this part from main programme to a seperate function
%%% Add b to select other wave shapes than only Sine
%--------------------------------------------------------------------------
%
% --------------------------------------
% |\   b   |                    |                   | 
% |   \     |       0.5         |    > 0.5        |
% |  r   \ |                     |                   |
% |-------------------------------------
% | 0.5    |       sine       |   acc. skewed|
% |         |                   |   (sawtooth)  |
% | ------------------------------------
% | >0.5  | vel. skewed   |  vel. +          |
% |         | (2nd Stokes) |   acc. skewed|
% |-------------------------------------


%% Calculation Tc  (as in 1C05 model)

if r==0.5  % Sine waves  or acceleration skewed (sawtooth)
       % Sine waves        
       if unet==0
           Tc=0.5*t;
       else
           if uwt<=unet*cos(Ang/180*pi)
               Tc=t;
           elseif uwc<=-unet*cos(Ang/180*pi)
               Tc=0;
           else
               Tc=t/pi*acos((-2*unet*cos(Ang/180*pi))/(uwc+uwt));
           end
       end   
elseif r >0.485   % also include near 0.5 measures
        % alternative way Tc and Tt; use sine function
        Tcw = t/pi*acos((-0.5*(uwc+uwt)+0.5*(9*uwc^2+9*uwt^2-14*uwc*uwt)^0.5)/...
                    (2*(uwc-uwt)));
        Ttw = t-Tcw;
        if unet * cos(Ang/180*pi) >= 0
            if uwt <= (unet*cos(Ang/180*pi)) 
                Tc=t;
            else
                Tc=t-Ttw*(1.0-2*asin(unet*cos(Ang/180*pi)/uwt)/pi);
            end
        else
            if uwc <= (-unet*cos(Ang/180*pi))
                Tc = 0;
            else
                Tc=Tcw*(1.0-2*asin(-unet*cos(Ang/180*pi)/uwc)/pi);
            end
       end
else
     % unknown value of r 
     disp(['Error wctp: Unknown combination r and b.  r=',num2str(r),'  b=', num2str(b)]);     
     Tc=NaN; 
end
  
Tt=t-Tc;

%% Calculation Tcu, Tcd, Ttu, Ttd (2nd order Stokes)

if abs(Tc)  < eps  % check on zero value
    Tcu = 0;
else
    Tcu = Tc*acos(2*b-1)/pi; %(see SANTOSS report J. Malarkey)
end
   
if abs(Tt)  < eps  % check on zero value
    Ttd = 0;
else
    Ttd = Tcu * Tt / Tc;
end
  
Tcd = Tc - Tcu;
Ttu =  Tt - Ttd;





