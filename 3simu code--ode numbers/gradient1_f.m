function df = gradient1_f(f,h)
df = zeros(1,length(f));
for j = 2:length(f) 
df(j) = (f(j) - f(j-1))/h; %后向差分
end
df(1) = df(2);%前向差分
end