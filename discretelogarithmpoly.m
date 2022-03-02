for ii=1:2
   
%d= randi([100 110])
eps=0.00000001;
d = 5;
%choose the power
p1=GenerateRandomPolynomial(5,1,10)
p2=GenerateRandomPolynomial(4,1,10)
q1=GenerateRandomPolynomial(5,1,10)
q2=GenerateRandomPolynomial(4,1,10)
%m=13
%choose matrices M and H

A=GenerateRandomMatrix(d,d,-1000,1000);
B=GenerateRandomMatrix(d,d,-1000,1000);
W=GenerateRandomMatrix(d,d,-1000,1000);


%ma=tp.generate_exponent(20)
%disp('Now generating an instance to which discrete log will be applied');
U=otimes(otimes(Applypolynomial(p1,A),W),Applypolynomial(p2,B))
V=otimes(otimes(Applypolynomial(q1,A),W),Applypolynomial(q2,B))
%U1=otimes(otimes(maxpower(A,m),W),maxpower(B,n));
%V1=otimes(otimes(maxpower(A,q),W),maxpower(B,r));
%disp('Now applying the discrete log');
[mlog,nlog]=discretelogstickel(U,W,A,B);
Ka=otimes(Applypolynomial(p1,A),otimes(V,Applypolynomial(q1,B)));
%Ka1=otimes(maxpower(A,m),otimes(V,maxpower(B,n)));
mattack=round(mlog);
nattack=round(nlog);
p1;
p2;
Uattack=otimes(otimes(powmaxplus(A,mattack),W),powmaxplus(B,nattack));
%Uattack1=otimes(otimes(maxpower(A,mattack),W),maxpower(B,nattack));
Kattack=otimes(powmaxplus(A,mattack),otimes(V,powmaxplus(B,nattack)));
%Kattack1=otimes(maxpower(A,mattack),otimes(V1,maxpower(B,nattack)));
if  (abs(Kattack-Ka)<=eps)
    k=k+1;
    disp('yes')
else
    disp('no')
    %break
end
%if  (abs(Kattack1-Ka1)<=eps)
 %   l=l+1;
  %  disp('yes2')
%else
 %   disp('no2')
    %break
%end
%D=Ka-Kattack
%[n,n]=size(D);

end
broken=k