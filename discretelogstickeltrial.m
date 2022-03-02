%trial discrete logarithm
clc;
clear all;
%dimension
k=0;
%l=0;
for ii=1:100
   
%d= randi([100 110])
eps=0.00000001;
d = 5;
%choose the power
m=randi([d^2 3*(d^2)]);
n=randi([d^2 3*(d^2)]);
q=randi([d^2 3*(d^2)]);
r=randi([d^2 3*(d^2)]);
%m=13
%choose matrices M and H

A=GenerateRandomMatrix(d,d,-1000,1000);
B=GenerateRandomMatrix(d,d,-1000,1000);
W=GenerateRandomMatrix(d,d,-1000,1000);


%ma=tp.generate_exponent(20)
%disp('Now generating an instance to which discrete log will be applied');
U=otimes(otimes(powmaxplus(A,m),W),powmaxplus(B,n))
V=otimes(otimes(powmaxplus(A,q),W),powmaxplus(B,r));
%U1=otimes(otimes(maxpower(A,m),W),maxpower(B,n));
%V1=otimes(otimes(maxpower(A,q),W),maxpower(B,r));
%disp('Now applying the discrete log');
[mlog,nlog]=discretelogstickel(U,W,A,B);
Ka=otimes(powmaxplus(A,m),otimes(V,powmaxplus(B,n)));
%Ka1=otimes(maxpower(A,m),otimes(V,maxpower(B,n)));
mattack=round(mlog)
nattack=round(nlog)
m
n
Uattack=otimes(otimes(powmaxplus(A,mattack),W),powmaxplus(B,nattack))
%Uattack1=otimes(otimes(maxpower(A,mattack),W),maxpower(B,nattack));
Kattack=otimes(powmaxplus(A,mattack),otimes(V,powmaxplus(B,nattack)));
%Kattack1=otimes(maxpower(A,mattack),otimes(V1,maxpower(B,nattack)));
%if  ((Kattack-Ka)<=eps)&& (mattack~=m)
if  ((Kattack-Ka)<=eps) 
    k=k+1;
    disp('yes')
   
else
    disp('no')
    %break
end
end
%if  (abs(Kattack1-Ka1)<=eps)
 %   l=l+1;
  %  disp('yes2')
%else
 %   disp('no2')
    %break
%end



