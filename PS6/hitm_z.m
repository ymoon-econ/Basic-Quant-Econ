function zdis = hitm_z(z,PPP,a1)
%HITM_Z simulates a finite state autoregressive time series using the 
%   pre-specified grid points and the associated transition matrix.
%
%                ydisc = hitm_z(z, PPP,a1)
%
%   where 
%
%       z - a column vector of the grid points,
%       PPP - the associated transition matrix,
%       a1 - a column vector of random numbers drawn from the standard 
%            uniform distribution on the interval(0,1), and
%       ydisc - the simulated time series. 
%
%   The length of the time series is given by the size of a1, i.e. if a1 is 
%   an Nx1 column vector of the random numbers, the size of ydisc is Nx1. 
%
%   Examples:
%   N=100; ydisc = hitm_z(z, PPP,rand(N,1))
%   or 
%   b=rand(100,1); ydisc = hitm_z(z, PPP,b);
%
%   Note that the method assumes that (i,j)-th element of PPP is the 
%   conditional probability Prob(state(t+1)=i|state(t)=j), i.e.
%
%                 PPP(i,j)=Prob(state(t+1)=i|state(t)=j)
%
%   where i and j denote the numbers of the states. Therefore 
%
%           PPP(1,j)+PPP(2,j)+ ... +PPP(n,j)=1, for all j.
%   See also HITM_S.  
%   Damba Lkhagvasuren, June 2005

N=size(a1,1);
znum=size(PPP,1);

A=tril(ones(znum))*PPP; A(znum,:)=2;

des=zeros(N,1);

ainit=randperm(znum);
des(1,1)=ainit(1,1);
destemp=des(1,1);

for c_ount=2:N;
        
          if a1(c_ount,1)<=A(1,destemp);
            des(c_ount,1)=1;
          end ;
          
          for i=1:znum-1;
              if A(i,destemp)<a1(c_ount,1)
                  if A(i+1,destemp)>=a1(c_ount,1);
                      des(c_ount,1)=i+1;
                  end
              end
          end
          destemp=des(c_ount,1);
end;

zdis=z(des);
