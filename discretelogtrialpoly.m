%trial discrete logarithm
clc;
clear all;
%dimension
%for i=1:100
%d= randi([100 110])
d = 5;
%choose the power
m=randi([((d-1)^2)+1 3*d^2])
n=randi([((d-1)^2)+1 2*d^2])
%choose matrices M and H
M1=GenerateRandomMatrix(d,d, -1000,1000);
H=GenerateRandomMatrix(d,d, -1000,1000)
%compute the message of Alice
%disp('Computing Alice message by definition:');
p1=GenerateRandomPolynomial(d,1,10)
p2=GenerateRandomPolynomial(d,1,10)
Ma=Applypolynomial(p1,M1)
Mb=Applypolynomial(p2,M1)
A=binexoticpower(M1,H,m);

%disp('Computing Alice message in a smart way:')
%A1=exoticpowerfirst(M,H,m);

%Resa=A-A1;
%if norm(Resa)>eps
%    error('The smart way of computing exotic power is not good?');
%end 

%compute the message of Bob
%disp('Computing Bob message by definition:');
B=binexoticpower(M1,H,n);

%disp('Computing Bob message in a smart way:')
%B1=exoticpowerfirst(M,H,n);

%Resb=B-B1;
%if norm(Resb)>eps
%    error('The smart way of computing exotic power is not good?');
%end 


%mth adjoint power of H
Hm=binadjointpower(H,m);

%nth adjoint power of H
Hn=binadjointpower(H,n);

%computing Alice's key
Ka=max(max(max(A,Applypolynomial(p1,B)),Hm), otimes(B,Hm));

%computing Bob's key
Kb=max(max(max(Applypolynomial(p2,A),B),Hn), otimes(A,Hn));

Reskakb=Ka-Kb

%H=[1 7 2 5;-1 -2 2 4;3 4 2 2;-5 -10 10 0]



%ma= 28
%nb= 31
%M=GenerateRandomMatrix(22,22,-1000,1000)


%I =maxeins(d)
%F=add(I,H)
%V=add(otimes(M,add(I,H)),H)

%disp('Now generating an instance to which discrete log will be applied');
%B=otimes(V,powmaxplus(H,n));
%A=otimes(V,powmaxplus(H,m))
%disp('Now applying the discrete log');
%nlog=discretelogavfk(B,V,H)
%mlog=discretelogavfk(A,V,H)
%if mlog==m && nlog==n
    %disp('success')
%else
   % disp('failed')
    %break
%end
%d;
%n;
%m

