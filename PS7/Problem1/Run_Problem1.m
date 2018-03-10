clc;
clear;

global a0 beta sigma r w T e

% Parameters
beta=1;
sigma=2;
r=0.03;
w=1;

% Initial Condition
a0=0;

% Generating Wage Efficiency Sequence (Quadratic on BLS Data)
e=zeros(1,77);
q0=1.16687413663967;
q1=0.76251324045698;
q2=-0.00621297952227995;

for i=1:77
    e(1,i) = q0 + q1*i + q2*i^2; 
end

% Forward Iteration
iter1=75;
iter2=500;
tol=0.000000001;

for T=1:iter1
    
    al = -10;
    ah = 10;
    
    %afinal_1 = fzero(@func1, al);
    
    % Use Bisection to find root of func1.m
    for i=1:iter2
    
        [al_T1,~,~] = func1(al);
        [ah_T1,~,~] = func1(ah);
    
        if ah_T1*al_T1>0
             disp('Guess again!')
            break
        else
            am=(ah+al)/2; 
            [am_T1,~,~]=func1(am); 
    
            if am_T1>0 
                ah=am;
            else
                al=am; 
            end
     
        end
    
        if abs(al-ah)<tol
            break
        else
    
        end
    end
end

T_star = T;
[afinal,a_1,c_1] = func1(ah);

% Generating Plots
age1 = linspace(16,90,75);
age = linspace(16,90,77);
e_1 = e(:,1:77);
a_1 = a_1(:,1:77);
c_1 = c_1(:,1:75);

figure
subplot(1,3,1);
plot(age,e_1);
title('Income Profile')
xlabel('Age')


subplot(1,3,2);
plot(age,a_1);
title('Assets');
xlabel('Age')


subplot(1,3,3);
plot(age1,c_1);
title('Consumption');
xlabel('Age')
