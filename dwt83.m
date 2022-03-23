function duv = dwt83(st,v)
% this function computes for constant alpha the 
% derivative of the wt83 parameterization
% Input
%      ST    standart deviation of orbiital velocity
%      V     mean alongshore current
% Output
%      DUV  derivative  

tmp1 = (v./st).^2;
tmp2 = sqrt(1.16^2 + tmp1);

duv = st.*(tmp2 + tmp1./tmp2);
