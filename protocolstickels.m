clc;
clear all;
k=0
for i =1:100

d = 5;
%choose the power
%m=randi([d^2 3*(d^2)]);
%n=randi([d^2 3*(d^2)]);
%q=randi([d^2 3*(d^2)]);
%r=randi([d^2 3*(d^2)]);
%m=13
%choose matrices M and H

%A=GenerateRandomMatrix(d,d,-1000,1000)
%B=GenerateRandomMatrix(d,d,-1000,1000)
%W=GenerateRandomMatrix(d,d,-1000,1000)
A = [595   432  -120  -959  -755;444  -901   395   389   793;-387  -990   357   211   474;558  -482   897   307  -233;81   186    27  -246  -829]
B =[ 607   476   749  -244  -775;310  -496   -74   917  -677;-270   479   103  -939   394;-105  -661  -999   500    67;-133  -172  -163  -508   935]
W = [-831   931   287   790  -250;-990  -901  -411  -424  -588;229  -989  -474   106  -714;141  -376   454   666   688;-450   820  -997  -699  -858]
m= 68;
n= 72;
q=53;
r=53;


%ma=tp.generate_exponent(20)
%disp('Now generating an instance to which discrete log will be applied');
U=otimes(otimes(powmaxplus(A,m),W),powmaxplus(B,n))
V=otimes(otimes(powmaxplus(A,q),W),powmaxplus(B,r))
%U1=otimes(otimes(maxpower(A,m),W),maxpower(B,n));
%V1=otimes(otimes(maxpower(A,q),W),maxpower(B,r));
%disp('Now applying the discrete log');
%[mlog,nlog]=discretelogstickeloldversionmain(U,W,A,B);
Ka=otimes(powmaxplus(A,m),otimes(V,powmaxplus(B,n)))
Kb=otimes(powmaxplus(A,q),otimes(U,powmaxplus(B,r)))
%Ka1=otimes(maxpower(A,m),otimes(V,maxpower(B,n)));
mattack=101
nattack=
m
n
q
r
Uattack=otimes(otimes(powmaxplus(A,mattack),W),powmaxplus(B,nattack));
%Uattack1=otimes(otimes(maxpower(A,mattack),W),maxpower(B,nattack));
Kattack=otimes(powmaxplus(A,mattack),otimes(V,powmaxplus(B,nattack)));
%Kattack1=otimes(maxpower(A,mattack),otimes(V1,maxpower(B,nattack)));
%if  ((Kattack-Ka)<=eps)&& (mattack~=m)
UU=U-Uattack
KK=Ka-Kb
if  (abs(Kattack-Ka)<=eps) 
    k=k+1;
    disp('yes')
    break
   
else
    disp('no')
    %break
end
end
k