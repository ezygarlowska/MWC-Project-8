function [Ux] = undertow(E,Er,c,h,rho,Hrms)
%Untertow computes the magnitude of the undertow
htrough=h-Hrms/2; 
M=E./c+2*Er./c;
Ux=M./(rho*htrough);
end

