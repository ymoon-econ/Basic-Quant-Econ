clc;
clear;

global a0 h0 alpha beta delta theta sigma chi r T gama

% Parameters
alpha=0.75;
beta=1;
delta=0.025;
theta=0.8;
sigma=2;
r=0.03;
chi=0.50;   % chi=1/nu
%w=1;
gama=0.8;

% Initial Conditions
a0=0;
h0=1;

% Forward Iteration
iter1=75;
iter2=1;
tol=0.000000001;

for T=1:iter1
    
    ch=4.9;
    nh=0.8985;
    
    cl=4.9;
    nl=0.8985;
    
    % Modify bisection to adjust for higher dimensions
%     for i=1:iter2    
%         
%         [ah_T1,~,~,~,~] = func2fwd(ch,nh);
%         [al_T1,~,~,~,~] = func2fwd(cl,nl);
%     
%         cm = (cl + ch)/2;
%         nm = (nl + nh)/2;
%          
%         [am_T1_1,~,~,~,~]=func2fwd(cm,nh);
%         [am_T1_2,~,~,~,~]=func2fwd(cm,nl);
%         [am_T1_3,~,~,~,~]=func2fwd(ch,nm);
%         [am_T1_4,~,~,~,~]=func2fwd(cl,nm);
%         [am_T1_5,~,~,~,~]=func2fwd(cm,nm);
         
%      
%         if am_T1_1>0
%            ch=cm;
%         else
%            cl=cm;
%         end
%         
%         if am_T1_2>0
%            ch=nm;
%         else
%            cl=cm;
%         end
%            
%         if abs(cl-ch)<tol
%             break
%         else
%             
%         end
    
%     end
    
end

[afinal,a,c,n,h] = func2fwd(cl,nl); 

% Generating Plots
age = linspace(16,90,77);
age2 = linspace(16,90,75);
a_2 = a(:,1:75);
a_1 = a(:,1:77);
c_1 = c(:,1:77);
n_1 = n(:,1:77);
h_1 = h(:,1:77);

figure
subplot(2,2,1);
plot(age2,a_2);
title('Assets');

subplot(2,2,2);
plot(age,c_1);
title('Consumption');

subplot(2,2,3);
plot(age,n_1);
title('Labor');

subplot(2,2,4);
plot(age,h_1);
title('Human Capital');