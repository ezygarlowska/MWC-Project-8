function [R,Beta]=velocity_skewness_asymmetry(u,t)
uc=max(u);
ut=abs(min(u));
dt=t(2)-t(1);
a=gradient(u,dt);
ac=max(a);
at=abs(min(a));
R=uc/(uc+ut);
Beta=ac/(ac+at);