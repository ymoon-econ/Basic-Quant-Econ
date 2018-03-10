function F = func1brute(x)

global a0 a76 beta sigma r et et1 age

a = x;
a_p = a(1:(age-1),1);
a_p = [a0; a_p];

a_f = a(2:age,1);
a_f = [a_f; a76];

F = ((1+r).*a + et1 - a_f) - ((beta*(1+r))^(1/sigma)).*((1+r).*a_p + et - a);

