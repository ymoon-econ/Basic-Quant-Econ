function [afinal,a,c,n,h] = func2LBND(c0,n0)

global a0 h0 alpha beta delta theta sigma chi r T gama

c = ones(1,T);      % Placeholder for consumption sequence
n = ones(1,T+1);    % Placeholder for labor sequence
a = ones(1,T+2);    % Placeholder for asset sequence
h = ones(1,T+2);    % Placeholder for human capital sequence

% Given initial conditions

a(1,1) = a0;    % a0 given
h(1,1) = h0;    % h0 given

c(1,1) = c0;    % c0 guess
n(1,1) = n0;    % n0 guess

for t=2:T+2
    
    % Finding c_t+1 (Equation 13)
    c(1,t) = ((beta*(1+r))^(1/sigma))*c(1,t-1);

    % Solving for n_t+1 using nonlinear equation solver
    x=c(1,t-1);
    y=c(1,t);
    v=h(1,t-1);
    w=h(1,t);
    z=n(1,t-1);
    
    % Finding h_t+1 (Equation (14))
    h(1,t) = (1-delta)*h(1,t-1) + theta*((1-n(1,t-1))^alpha);
    
    % Finding a_t+1 (Equation (15))
    a(1,t) = (1+r)*a(1,t-1) + h(1,t-1)*n(1,t-1) - c(1,t-1);
    
  
    % Equation (16) 
    f = @(n) beta*(y^(-sigma))*n...
                + beta*(1-delta)*((1/(alpha*theta))*((1-n)^(1-alpha))*(((1-n)^(-gama)) + (y^(-sigma))*w))...
                - (1/(alpha*theta))*((1-z)^(1-alpha))*(((1-z)^(-gama)) + (x^(-sigma))*v);
    x0=0.5; % guess to solve nonlinear equation (16)
    ntp1 = fsolve(f,x0);       

%     % Alternate Equation (16)
%     f = @(n) beta*(x^(-sigma))*n...
%                 + beta*(1-delta)*((1/(alpha*theta))*((1-n)^(1-alpha))*((y^(-sigma))*w))...
%                 - (1/(alpha*theta))*((1-z)^(1-alpha))*((x^(-sigma))*v); 
    %ntp1 = real(ntp1);
    %ntp1 = fsolve(f,x0);
    
    % Checking corner solution
        if ntp1>0    
             n(1,t)=ntp1;
        else
             n(1,t)=0;
        end
      
end

afinal = a(1,T+2);

