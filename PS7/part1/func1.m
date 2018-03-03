function [afinal, a, c] = func1(a1)            % a1 is what I guess


global a0 beta sigma r w e T

a = zeros(1, T+2);
a(1,1) = a0;
a(1,2) = a1;

c = zeros(1, T);


for t=3:T+2
    
    a(1,t) = (1+r)*a(1,t-1) + w*e(1,t-1)...
       - ((beta*(1+r))^(1/sigma))*((1+r)*a(1,t-2) + w*e(1,t-2) - a(1,t-1));
    
    c(1,t-2) = (1+r)*a(1,t-2) + w*e(1,t-2) - a(1,t-1);    

end

afinal = a(1,T+2);
