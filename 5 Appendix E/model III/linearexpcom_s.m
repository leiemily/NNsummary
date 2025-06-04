function ss = linearexpcom_s(ak,TT,tt,step_maxvalue,per)
  ss = zeros(1,length(tt));
  for k = 1:length(ak)
      ss = ss + ak(k)*(1-cos(2*pi*tt/TT(k)));
  end
%   ss = exp(ss - (abs(log(step_maxvalue)) + step_maxvalue));
  ss = exp(ss - (max(ss)-abs(log((1-per)*step_maxvalue))));
%   ss = exp(ss - (step_maxvalue-abs(log(step_maxvalue))));
end