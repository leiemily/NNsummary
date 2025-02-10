function ss = linearexpcom_s(ak,TT,tt,step_maxvalue)
  ss = zeros(1,length(tt));
  for k = 1:length(ak)
      ss = ss + ak(k)*(1-cos(2*pi*tt/TT(k)));
  end
%   ss = exp(ss - (abs(log(step_maxvalue)) + step_maxvalue));
  ss = exp(ss - (abs(log(step_maxvalue)) + max(ss)));
end