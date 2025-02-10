function df2 = gradient2_f(f,h)
df2 = zeros(1,length(f));
for j = 2:(length(f)-1) 
df2(j) = (f(j-1) - 2*f(j) + f(j+1))/h^2; %三点中心差分
end
df2(1) = (f(1) - 2*f(2) + f(3))/h^2;%三点前向差分
df2(end) = (f(end-2) - 2*f(end-1) + f(end))/h^2;%三点后向差分
end

% df2 = zeros(1,length(f));
% for j = 4:length(f)
% df2(j) = (2*f(j) - 5*f(j-1) + 4*f(j-2) - f(j-3))/h^2; %四点向后差分
% end
% df2(1) = (f(1) - 2*f(2) + f(3))/h^2;%三点前向差分
% df2(2) = (f(1) - 2*f(2) + f(3))/h^2; %三点中心差分
% df2(3) = (f(2) - 2*f(3) + f(4))/h^2; %三点中心差分
% end