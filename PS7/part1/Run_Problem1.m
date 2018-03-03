clc;
clear;

%Parameters
global a0 beta sigma r w T e

a0=0;
beta=1;
sigma=2;
r=0.03;
w=1;
e=1;

e=zeros(1,77);
q0=1.16687413663967;
q1=0.76251324045698;
q2=-0.00621297952227995;

for i=1:77
    e(1,i) = q0 + q1*i + q2*i^2;
    
end


age = linspace(16,90,75);
age1 = linspace(16,92,77);

%plot(age,e)

iter1=75;
tol=0.000000001;
iter2=500;

for T=1:iter1
    
    al = -10;
    ah = 10;
    
    for i=1:iter2
    
        [al_T1,~,~] = func1(al);
    
        [ah_T1,~,~] = func1(ah);
    
        if ah_T1*al_T1>0
             disp('Bad guess, recheck your guess')
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

a_1 = a_1(:,1:75);
plot(age,a_1)
%plot(age,c_1)