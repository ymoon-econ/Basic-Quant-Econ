% written by Priyam

global sigmae rho m alpha delta eta b beta

sigmae=0.0034;
rho=0.98;
m=2;                    % arbitrary choice between 2 and 4
alpha=0.72;
delta=0.1;
eta=0.72;
b=0.4;
beta=0.99;
n=30;                   % Number possible realizations of the process of z
%Discretizing an AR(1) process of productivity shocks
%logZ(t)=rho*logZ(t-1)+sigmae*Et 
%where Et is the error term following a normal distribution (0,1)
% Then, sigmae*Et follows a normal distribution (0, sigmae^2)


%Tauchen Method
zn=m*(sigmae)./sqrt(1-rho^2);
z1=-m*(sigmae)./sqrt(1-rho^2);
Z=linspace(z1,zn,n);
d=Z(2)-Z(1);                    % Z has equal distance
ptau=ones(n,n);                 % Transition Matrix
for k=2:n-1
    for j=1:n
       ptau(j,k)=normcdf((Z(k)+d/2-rho*Z(j))/sigmae)-normcdf((Z(k)-d/2-rho*Z(j))/sigmae);
    end
end

for j=1:n
    ptau(j,1)=normcdf((Z(1)+d/2-rho*Z(j))/sigmae);
    ptau(j,n)=1-normcdf((Z(1)-d/2-rho*Z(j))/sigmae);
end
%Rowenharst Method
p=(1+rho)/2;
q=p;
zhi=sqrt(n-1)*((sigmae)./sqrt(1-rho^2));
P_Rouw=[ p  (1-p);
        (1-q) q];
    
    for i_R=2:n-1
    a1R=[P_Rouw zeros(i_R, 1); zeros(1, i_R+1)];
    a2R=[zeros(i_R, 1) P_Rouw; zeros(1, i_R+1)];
    a3R=[zeros(1,i_R+1); P_Rouw zeros(i_R,1)];
    a4R=[zeros(1,i_R+1); zeros(i_R,1) P_Rouw];
    P_Rouw=p*a1R+(1-p)*a2R+(1-q)*a3R+q*a4R;
    P_Rouw(2:i_R, :) = P_Rouw(2:i_R, :)/2;
    end
    
P_Rouw=P_Rouw';

for i_R = 1:n
    P_Rouw(:,i_R) = P_Rouw(:,i_R)/sum(P_Rouw(:,i_R));
end
%%%%%

%Generate data from Rowenharst method
a1=randn(n);
Zsim=fsolve(hitm_z,Z,P_Rouw,a1);
