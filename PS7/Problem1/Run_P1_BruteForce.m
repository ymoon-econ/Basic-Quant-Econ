clc;
clear;

% Parameters
beta = 1;
sigma = 2;
r = 0.03;
age = 75;
%w = 1;

% Initial and Boundary Conditions
a0 = 0;
a76 = 0;

% Generating Wage Efficiency Sequence (Quadratic on BLS Data)
e = zeros(1,age);

q0 = 1.16687413663967;
q1 = 0.76251324045698;
q2 = -0.00621297952227995;

for i=1:age
    e(1,i) = q0 + q1*i + q2*i^2; 
end

e = e';
et = e(1:age,1);
et1 = e(2:age,1);
et1 = [et1; 0];

% Solve System of Euler Equations
x0 = ones(age,1);
a_sol = fsolve(@func1brute, x0);

% Generate Sequence of Assets and Consumptions
a = [a0; a_sol];

ac = a(1:age,:);
c = (1+r).*ac + e - a_sol;

% Generate Plots
years0 = linspace(16,90, age);
years = linspace(16,90,(age+1));

figure
subplot(2,2,1)
plot(years0,e);
title('Income Profile')

subplot(2,2,2);
plot(years,a);
title('Assets');

subplot(2,2,3);
plot(years0,c);
title('Consumption');
 