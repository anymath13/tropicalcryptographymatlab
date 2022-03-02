%trial discrete logarithm
clc;
clear all;
%dimension
k=0;
%l=0;
for ii=1:100
   
%d= randi([100 110])
eps=0.000001;
d = 50;
%choose the power
m=randi([d^2 3*(d^2)]);
n=randi([d^2 3*(d^2)]);
q=randi([d^2 3*(d^2)]);
r=randi([d^2 3*(d^2)]);
%m=13
%q=53
%r=53
%choose matrices M and H

A=GenerateRandomMatrix(d,d,-10^5,0);
B=GenerateRandomMatrix(d,d,-10^5,0);
W=GenerateRandomMatrix(d,d,0,10^5);
%A = [595   432  -120  -959  -755;444  -901   395   389   793;-387  -990   357   211   474;558  -482   897   307  -233;81   186    27  -246  -829]
%B =[ 607   476   749  -244  -775;310  -496   -74   917  -677;-270   479   103  -939   394;-105  -661  -999   500    67;-133  -172  -163  -508   935]

%W = [-831   931   287   790  -250;-990  -901  -411  -424  -588;229  -989  -474   106  -714;141  -376   454   666   688;-450   820  -997  -699  -858]

%m= 68;
%n= 72;



%ma=tp.generate_exponent(20)
%disp('Now generating an instance to which discrete log will be applied');
U=otimes(otimes(maxpower(A,m),W),maxpower(B,n));
V=otimes(otimes(maxpower(A,q),W),maxpower(B,r));
%U1=otimes(otimes(maxpower(A,m),W),maxpower(B,n));
%V1=otimes(otimes(maxpower(A,q),W),maxpower(B,r));
%disp('Now applying the discrete log');
[mlog,nlog]=discretelogstickeloldversionmain(U,W,A,B);
[qlog,rlog]=discretelogstickeloldversionmain(V,W,A,B);
Ka=otimes(maxpower(A,m),otimes(V,maxpower(B,n)));
%Ka1=otimes(maxpower(A,m),otimes(V,maxpower(B,n)));
mattack=round(mlog);
nattack=round(nlog);
%mattack=35
%nattack=93
%qattack=randi(100)
%rattack=randi(100)
qattack=round(qlog);
rattack=round(rlog);
m;
n;
q;
r;
%Uattack=otimes(otimes(maxpower(A,mattack),W),maxpower(B,nattack));
%Vattack=otimes(otimes(maxpower(A,qattack),W),maxpower(B,rattack));
Kattack=otimes(maxpower(A,mattack),otimes(V,maxpower(B,nattack)));
%Kattack1=otimes(maxpower(A,qattack),otimes(U,maxpower(B,rattack)));
%if  ((Kattack-Ka)<=eps)&& (mattack~=m)
%UU=U-Uattack;
%VV=V-Vattack;
if  (abs(Kattack-Ka)<=eps) 
    k=k+1;
    disp('yes')
    %break;
   
else
    disp('no')
    %break
end
end
k

%if  (abs(Kattack1-Ka1)<=eps)
 %   l=l+1;
  %  disp('yes2')
%else
 %   disp('no2')
    %break
%end



