%Kim Ruhl May 20, 2004, this version: April 2008
%(last modification: modified to accept natural units for working age pop)
%This file loads the data and parameters for the model and calls the
%SolveModel function. To run this script, make sure this file is in your
%working directory, or the directory is in the MATLAB path. Then simply
%type "depressions" without quotes, to run the script.

%The file "paramBase.txt" is a column vector with the following ordering of parameters:
%beta, gamma, delta, alpha, growth factor for TFP on BGP, growth factor for populatin on BGP
%and the initial capital stock.

% The file "dataBase.txt" is a Tx6 matrix where rows are years, and columns are:
% TFP, working age pop., available hours, consumption tax rate, labor tax rate, capital tax rate.

%Clear the workspace and load parameters and data.  These files can be
%created by Excel using the "save as tab delimited text" option.
clear;
param=load('paramBase.txt');
data=load('dataBase.txt');


%delcare the global parameters
global BETA;    global GAMMA;   global DELTA;   global ALPHA;   global G;
global K0;      global T;       global ETA;

%assign them values from the file
BETA = param(1);    %discount factor    
GAMMA = param(2);   %utility function parameter
DELTA = param(3);   %depreciation rate
ALPHA = param(4);   %capital's share in production
G = param(5);       %TFP growth factor in balanced growth path
ETA = param(6);     %population growth factor in balanced growth path 
K0 = param(7);      %initial capital stock

T= size(data,1);    %number of periods
%assign the values from the file
a = data(:,1);
pop = data(:,2);
lbar = data(:,3);
tauc = data(:,4);
taul = data(:,5);
tauk = data(:,6);

%create TFP, population, and available hours that grow at constant rates
%the solution the model with constant series should be easy to solve
acst(1,1)=a(1,1);
popcst(1,1)=pop(1,1);
lbarcst(1,1)=lbar(1,1);
for i=2:length(a)
    acst(i,1)   =acst(i-1,1)*G;
    popcst(i,1) =popcst(i-1,1)*ETA;
    lbarcst(i,1)=lbarcst(i-1,1)*ETA;
end

%start with no taxes
tauccst = zeros(T,1);
taulcst = zeros(T,1);
taukcst = zeros(T,1);

%Use initial capital stock, and guess that capital grows at the BGP rate
x(1,1) = G^(1/(1-ALPHA))*ETA*K0;
for i = 2:T-1
    x(i,1) = G^(1/(1-ALPHA))*ETA*x(i-1);
end

%guess that leisure is (1-GAMMA)*(1-taul(i)) of available time
x(T:2*T-1,1) = (1-GAMMA).*(1.-taulcst).*lbarcst;

%taking logs guarantees the variables will be positive
x=log(x);

%Call SolveModel, and check for success.  Start with a model with constant growth rates (lam=1) 
%and move to the version you want to solve (lam=0). Use the
%solution from the model just solved as a begining guess for the next
%model.  If the routine is not finding solutions, smaller lam steps may help.

for lam = 1.0:-0.2:0.2
    fprintf('=========Lambda = %3.2f ==============\n', lam);
    [err,x] = solveModel( 0, x, lam*acst+(1-lam)*a, lam*popcst+(1-lam)*pop, lam*lbarcst+(1-lam)*lbar,...
                      lam*tauccst+(1-lam)*tauc, lam*taulcst+(1-lam)*taul, lam*taukcst+(1-lam)*tauk);
    if err == 0  
        fprintf('Model solved successfully.\n');
    else
        fprintf('Could not find model solution.\n');
    end
end

%ready to solve the version of the model that uses the data
fprintf('=========Lambda = 0.0 ==============\n');
[err,x] = solveModel(1,x,a,pop,lbar,tauc,taul,tauk);



        
