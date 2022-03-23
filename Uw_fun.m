function Uw = Uw_fun(h,Hrms,T)

k = k_fun(T,h);
Uw = (pi*Hrms)./(T.*sinh(k.*h));

end 