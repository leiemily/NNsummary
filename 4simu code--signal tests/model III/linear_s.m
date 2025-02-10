function [ss,ss_slope,ss_steplen] = linear_s(dt,discreteT_length,step_num,step_maxvalue,step_minlen,slope_max,slope_min)
t = 0:dt:(discreteT_length-1)*dt;
r1 = step_minlen/dt;
r2 = discreteT_length - r1;
ss = zeros(1,discreteT_length);
ss_slope = (slope_max+1)*ones(1,step_num-1);
ss_steplen = zeros(1,2*step_num-1);
while ~isempty(find(ss_steplen < step_minlen, 1)) || ~isempty(find(abs(ss_slope) > slope_max, 1)) || ~isempty(find(abs(ss_slope) < slope_min, 1))
    rng("shuffle");
    ss_r = sort(randi([r1,r2],1,2*(step_num-1)));
    ss_r = [1,ss_r,discreteT_length];
    step_value = unifrnd(0,step_maxvalue,1,step_num);
    for i = 1:step_num-1
        ss(ss_r(2*i)) = step_value(i);
        ss(ss_r(2*i+1)) = step_value(i+1);
        ss_steplen(2*i-1) = (ss_r(2*i)-ss_r(2*i-1))*dt;
        ss_steplen(2*i) = (ss_r(2*i+1)-ss_r(2*i))*dt;
        ss_slope(i) = (step_value(i+1)-step_value(i))/(ss_r(2*i+1)-ss_r(2*i))*dt^(-1);
    end
    ss_steplen(end) = (ss_r(end)-ss_r(end-1))*dt;
    ss(ss_r(1)) = step_value(1);
    ss(ss_r(end)) = step_value(end); 
end
ss = interp1(t(ss_r),ss(ss_r),t);
end