%protocol trial
%running protocol
eps=0.000000001;
%for i=1:100
%d= randi([40 50]) 
%d=randi([100 150])
%d=randi([2 30]);
d=3;
%choose an integer used by Alice
m=randi([((d-1)^2)+1 d^2]);
%m=[2^10 2^200]
m=6;
n=7;
%m= 449
%choose an integer used by Bob
n=randi([((d-1)^2)+1 d^2]);

%choose matrices M and H
%M=GenerateRandomMatrix(d,d,-1000,1000);
M= [[526, -1, 817]; [323, -757, -475]; [244, -9, -731]]
H= [[-971, 577, 15]; [-558, -816, -272]; [-734, 106, -413]]
%H=generatespecialmatrix(d);
%=[8 7 2;10 3 6;-10 -1 3]
%H=[0 -3 -5;-1 -2 2;1 -3 -4]
%H=GenerateRandomMatrix(d,d,-1000,1000);

%H=maxplussprand(d,d,10)
%H=[1 7 2 5;-1 -2 2 4;3 4 2 2;-5 -10 10 0];


%compute the message of Alice
disp('Computing Alice message by definition:');
A=newexoticpower(M,H,m);

disp('Computing Alice message in a smart way:')
%A1=exoticpowerfirst(M,H,m);

%Resa=A-A1;
%if norm(Resa)>eps
%    error('The smart way of computing exotic power is not good?');
%end 

%compute the message of Bob
disp('Computing Bob message by definition:');
B=newexoticpower(M,H,n);

%disp('Computing Bob message in a smart way:')
%B1=exoticpowerfirst(M,H,n);

%Resb=B-B1;
%if norm(Resb)>eps
%    error('The smart way of computing exotic power is not good?');
%end 

disp('computing H^m');
%mth adjoint power of H
Hm=newadjointpower(H,m);
disp('computing H^n');
%nth adjoint power of H
Hn=newadjointpower(H,n);

disp('computing Alice key');
Ka=max(max(max(A,B),Hm), otimes(B,Hm));

disp('computing Bob key');
Kb=max(max(max(A,B),Hn), otimes(A,Hn));