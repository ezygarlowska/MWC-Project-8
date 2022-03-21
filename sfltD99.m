%-------------------------------------------------------------------------
%   COMPUTATION OF THE SHEET-FLOW LAYER THICKNESS USING THE FORMULA OF
%   DOHMEN-JANSSEN (1999)
%
%   Created by Jebbe van der Werf on 26-09-2006
%   Modified by Jebbe van der Werf on 28-06-2007  
%   20 oct 2009, version 2.04, Dominic van der A/René Buijsrogge
%   Correction sheet flow layer thickness for fine sand due to increased
%   fine sand mobile roughness:
%   Now interpolation between D50=0.13 and 0.20mm (was D50=0.13 and 0.21mm)
%   Now SFLT = D50*25.08*Sw (was SFLT=D50*35*Sw)
%   Copyright 2009 University of Twente, Water Engineering and Management, Enschede, The Netherlands
%-------------------------------------------------------------------------

function [SFLTc,SFLTt]=sfltD99(D50,Sc,St,Swc,Swt,sett)
fs_const=25.08; %fine sand constant
if strcmp(sett(3),'W')==1
    %using the wave-alone bed shear stresses
    for i=1:length(D50)
        if D50(i)<=0.00015
            SFLTc(i,:)=D50(i)*fs_const*Swc(i);
            SFLTt(i,:)=D50(i)*fs_const*Swt(i);    
        elseif D50(i)>=0.00020
            SFLTc(i,:)=D50(i)*13*Swc(i);
            SFLTt(i,:)=D50(i)*13*Swt(i);
        else
            %linear interpolation between values for D50=0.15 mm and D50=0.20
            %mm
            SFLTc(i,:)=D50(i)*Swc(i)*(fs_const+(D50(i)-0.00015)*(13-fs_const)/(0.00020-0.00015));
            SFLTt(i,:)=D50(i)*Swt(i)*(fs_const+(D50(i)-0.00015)*(13-fs_const)/(0.00020-0.00015));
        end
    end
elseif strcmp(sett(3),'WC')==1
    %using the total bed shear stresses
    for i=1:length(D50)
        if D50(i)<=0.00015
            SFLTc(i,:)=D50(i)*fs_const*Sc(i);
            SFLTt(i,:)=D50(i)*fs_const*St(i);    
        elseif D50(i)>=0.00020
            SFLTc(i,:)=D50(i)*13*Sc(i);
            SFLTt(i,:)=D50(i)*13*St(i);
        else
            %linear interpolation between values for D50=0.15 mm and D50=0.20
            %mm
            SFLTc(i,:)=D50(i)*Sc(i)*(fs_const+(D50(i)-0.00015)*(13-fs_const)/(0.00020-0.00015));
            SFLTt(i,:)=D50(i)*St(i)*(fs_const+(D50(i)-0.00015)*(13-fs_const)/(0.00020-0.00015));
        end
    end
end