function ss = linearcom_s(ak,TT,tt,Numberof_Tp)
  ss = zeros(1,length(tt));
  for k = 1:length(ak)
      ss = ss + 1/Numberof_Tp*ak(k)*(1-cos(2*pi*tt/TT(k)));
  end
end