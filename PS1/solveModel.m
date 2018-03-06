%Kim Ruhl, May 20, 2004 this version: April 2008 
%(last modification: modified to accept natural units for working age pop)
%This function solves the system of equations specified in the subfunction EvalSystem.  
%There are T intratemperal FOCs and T-1 intertemperal FOCs, for a 2T-1 system of equations.

%The function uses Newton's method to solve the system of equations.
%The results are written to the file "output.xls"
%This function is called by depressions.m

function [err, sol] = solveModel(flag, x, a, pop, lbar, tauc, taul, tauk)
err = 0;
global BETA;    global GAMMA;   global DELTA;   global ALPHA;   global G;
global ETA;     global K0;      global T;       

ITR_MAX = 250;              %Maximum interations for Newton's method
EPSILON = 1.0e-10;          %Stopping value for Newton's method
J_STEP = 1.0e-08;           %Step size for computing numerical derivatives

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If the model does not solve, you can try changing these values.          
N_ITER = 3;                 %Number of iterations before N_STEP = 1.0     
N_STEP = 0.5;               %Begining step size in Newton's Method        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fx = zeros(2*T-1,1);            %value of the equations (what we want to be zero)
j = zeros(2*T-1,2*T-1);         %jacobian matrix
j_inv = zeros(2*T-1,2*T-1);     %inverse of the jacobian matrix

iteration = 0;
loss_fx = 1;
fx = EvalSystem(x, a, pop, lbar, tauc, taul, tauk);

%Newton's method
while (loss_fx >= EPSILON) & (iteration <= ITR_MAX)     %check: are we close enough to zero? have we tried too many times?
    j = Jacobian(fx, x, J_STEP, a, pop, lbar, tauc, taul, tauk);
    j_inv = inv(j);
    if( iteration >= N_ITER) N_STEP = 1.0; end              %Once things are going well, start taking bigger steps
    x = x - N_STEP.*(j_inv*fx);                             %update guess     
    fx = EvalSystem(x, a, pop, lbar, tauc, taul, tauk);
    
    iteration = iteration + 1;
    loss_fx = max(abs(fx));                            %how big is the largest error?
    fprintf('iteration %5.2f \t',iteration);
    fprintf('error %10.6e \n', loss_fx);   
end

%check for failure
if iteration > ITR_MAX
    fprintf('Max number of iterations exceeded.');
    err = 1;
end
sol = x;      %return the solution 

%if flag =  1 then compute results and output them to a file
if flag == 1
    %since we logged everything, now we turn it back
    xk(1) = K0;
    xk(2:T,1) = exp(x(1:T-1));
    xk(T+1) = xk(T)*G^(1/(1-ALPHA))*ETA;  
    xl = exp(x(T:2*T-1,1));
    
    %Using the correct capital and leisure values, find output per worker,
    %investment/GDP, hours worked per week, consumption/GDP, K/Y and the
    %interest rate.  This section can be modified to output any variables
    %of interest.
    
    results = zeros(T,6);
    y = a.*xk(1:T,1).^ALPHA .* (lbar - exp(x(T:2*T-1,1))).^(1-ALPHA); %GDP at factor prices
    c = y - xk(2:T+1,1) + (1-DELTA).*xk(1:T,1);
    y = y + c.*tauc;  %GDP at market prices
    
    %results(:,1) = y./(pop.*1000);
     results(:,1) = y./(pop);
    results(:,2) = (xk(2:T+1,1)-(1-DELTA).*xk(1:T,1))./y;
    results(:,3) = (lbar - xl)./lbar;
    results(:,4) = (y-xk(2:T+1,1)+(1-DELTA)*xk(1:T,1))./y;
    results(:,5) = xk(1:T,1)./y;
    results(:,6) = ALPHA.*a.*xk(1:T,1).^(ALPHA-1).*(lbar-xl).^(1-ALPHA)-DELTA;

    %write the data out to a file
    save('output.xls','results','-ascii', '-tabs');
    fprintf('Progam finished. Data written to file. \n');
end

%--------------------------------------------------------------------------
%----------------------------------EvalSystem------------------------------
%--------------------------------------------------------------------------
function f = EvalSystem(x, a, pop, lbar, tauc, taul, tauk)
%This function evaluates the system of equations.  It is called by Jacobian
%and SolveModel

%use the global parameter values
global BETA;    global GAMMA;   global DELTA;   global ALPHA;   global G;
global ETA;     global K0;      global T;       

f = zeros(size(x,1),1);

xk = zeros(T+1,1);      %capital stocks
xl = zeros(T,1);        %leisure 
n = zeros(T,1);         %labor
y = zeros(T,1);         %output
c = zeros(T,1);         %consumption
w = zeros(T,1);         %wages

%convert data from logs
xk(1) = K0;
xk(2:T,1) = exp(x(1:T-1));
xk(T+1) = G^(1/(1-ALPHA))*ETA*xk(T);

xl = exp(x(T:2*T-1,1));
%compute labor
n = lbar - xl;

%compute output, consumption and wages
y = a.*xk(1:T,1).^ALPHA.*n.^(1-ALPHA);
c = y - xk(2:T+1,1) + (1-DELTA).*xk(1:T,1);
w = (1-ALPHA).*a.*xk(1:T,1).^ALPHA.*n.^(-ALPHA);

%Compute the intratemperal FOC. The division by y scales the equation. 
f(1:T,1) = (w.*xl.*(1-taul) - (1-GAMMA)/GAMMA.*c.*(1+tauc))./y;

%Compute the intertemporal FOC.
f(T+1:2*T-1,1) = c(2:T,1)./c(1:T-1,1).*(1+tauc(2:T,1))./(1+tauc(1:T-1,1)) - ...
                BETA.*(1+(1-tauk(2:T,1)).*(ALPHA*a(2:T,1).*xk(2:T,1).^(ALPHA-1).*n(2:T,1).^(1-ALPHA) - DELTA));

%--------------------------------------------------------------------------
%----------------------------------Jacobian--------------------------------
%--------------------------------------------------------------------------
function jac = Jacobian(fx, x, step, a, pop, lbar, tauc, taul, tauk)
%This function computes the Jacobian matrix of the function in EvalSystem at x. 
%The argument step is the size of the step used in numerical differentiation, 
%the argument fx is the value of the system evaluated at x.  The Jacobian
%is returned in the matrix jac.

%Create matrix of correct size to hold function values
n = size(x,1);
fstep = zeros(n,1);

%For each function variable, increment the variable by step, recompute
%the function value, f(x+step), and record [f(x+step)-f(x)]/step.  
for c = 1:n
    x(c) = x(c) + step;
    fstep = EvalSystem(x, a, pop, lbar, tauc, taul, tauk);
    jac(:,c) = (fstep-fx)./step;  
    x(c) = x(c) - step;
end

