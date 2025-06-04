function df = gradient1_s(f,h)
df = zeros(1,length(f));
for j = 2:length(f)-1 
df(j) = (f(j+1) - f(j-1))/(2*h); %中心差分
end
df(1) = (f(2) - f(1))/h;%前向差分
df(end) = (f(end) - f(end-1))/h;%后向差分
end