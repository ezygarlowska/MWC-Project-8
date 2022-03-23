function [Ur] = Ursell(k,h,Hrms)
% This function computes the Ursell number from k, h and Hrms
    a=0.5*sqrt(2).*Hrms;
    Ur=(3*a.*k)./(4*(k.*h).^3);
end
