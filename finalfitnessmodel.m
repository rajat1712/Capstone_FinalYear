function F= finalfitnessmodel(in_para, alpha)
s = fotf('s')  % FOPID transfer function
kp=in_para(1,1);ki=in_para(1,2);kd=in_para(1,3);lambda=in_para(1,4); mu=in_para(1,5);% Getting corresponding values from matrix
sys= 761.44735/(s^2.0971+4.0972*s^1.0036+3.8842);  % Initial system
hs= ((s^(lambda)*kp)+ki+(s^(mu+lambda))*kd)/(s^(lambda));  % FOPID General Equation
result= sys*hs/(1+sys*hs);
[y,t]=step(result);
S = stepinfo(y,t); % Getting all parameters like rising time, maximum overshoot, delay time, peak time etc.
error = sum(abs(y-1)); % Subtracting 1 as we're taking unit step input and summing because of integration
F = alpha*error + (1-alpha)*(S.RiseTime+S.SettlingTime);  % Final value //u(t) missing//
