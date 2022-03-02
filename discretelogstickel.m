function [m,n]=discretelogstickel(U,W,A,B)
% Finds k such that A=V F^k
% works when lambda(F) is not equal to 0
% currently only for finite entries

%eps=10^{-9}
eps=0.000000001;

%disp('Discretelog starts here.');
%A
%V
%F

[d1,d2]=size(A);




%We first check if any of the first (n-1)^2 powers work.
%for higher powers, the method is based on the ultimate periodicity 
%of critical columns

%disp('compute lambda(A) and lambda(B)');

[chi1,x1,criticalcycle1]=policyIteration(A);
lambda1=max(chi1);
[chi2,x2,criticalcycle2]=policyIteration(B);
lambda2=max(chi2);
 if abs(lambda1)<eps
     error('The mcm of the third argument is equal to 0. The method does not work in this case.');
 end


AA=A-lambda1;
BB=B-lambda2;

e1=eigvector(AA);
e2=eigenvector(BB);

%disp('Now finding a critical cycle in A');
cycle1 =critical(AA,e1);
[~,l1]=size(cycle1);
c1=l1-1;
zz1=zeros(d1,1);
for i=1:c1
    zz1(cycle1(i))=1;
end
%zz
%disp('Now finding a critical cycle in B');
cycle2 =critical(BB,e2);
[~,l2]=size(cycle2);
c2=l2-1;
zz2=zeros(d1,1);
for i=1:c2
    zz2(cycle2(i))=1;
end
%disp('Now using the CSR to catch the power');


U1=maxfloyd(maxpower(AA,c1));
R=row(U1,zz1);


S1=-inf*ones(d1,d1);
for i=1:c1   
S1(cycle1(i),cycle1(i+1))=AA(cycle1(i),cycle1(i+1));
end
S1;



U2=maxfloyd(maxpower(BB,c2));
C=column(U2,zz2);


S2=-inf*ones(d1,d1);
for i=1:c2   
S2(cycle2(i),cycle2(i+1))=BB(cycle2(i),cycle2(i+1));
end
S2;

for j=1:c2
Y2(:,j)=U(:,cycle2(j));
end
Y2;

newflag=0; 
 for t1=1:c1 
     for t2=1:c2
 t1;
 t2;
 W1=otimes(binmaxpower(S1,t1),R);
 W2=otimes(C,binmaxpower(S2,t2));
 Ucand=otimes(W1,otimes(W,W2)); %+lambda1 + lambda2;
 %Y
 
 
 
 X=-inf*ones(c1,c2);
 for i=1:c1
     for j=1:c2 
         X(i,j)=Ucand(cycle1(i),cycle2(j));
     end
 end
 
 X;
 
  Y=-inf*ones(c1,c2);
 for i=1:c1
     for j=1:c2 
         Y(i,j)=U(cycle1(i),cycle2(j));
     end
 end
 
 Y;
 
 
 Z=Y-X;
 % will not be good if we have -inf entries..
 Mu1=Z(1,1);
 
 Z=Z-Mu1;
 if norm(Z)<eps
     newflag=1;
 break
 end

     end
     % we need to find m,n such that
     %m*l2*s1 + n*l1*s2=(Mu1*l1*l2)
     % we need to find a pair (m,n)
     %lambda1=s1/l1, lambda2=s2/l2
     % m in the range from 0 to floor((Mu*l1)/s1) for every m
     % in this range, check if (Mu1*l1*l2-m*l2*s1)/(l1*s2) is an integer
%  then it is n

    %power=-1;
    %power =k
%else
s1=lambda1*c1;
s2=lambda2*c2;



flag=0;
for m=d1^2:floor((Mu1*c1)/s1)
   c1;
   c2;
   t1;
   t2;
    m;
    m1=(m-t1)/c1;
    if (abs(m1-round(m1))<=eps)
        %disp('m is good, now checking n');
        n=(Mu1*c1*c2-m*c2*s1)/(c1*s2);
        n1=round(n);
        n2=(n1-t2)/c2;
        if (abs(n2-round(n2))<=eps)
           if (abs(n-n1)<= eps)
             % disp('n is also good, the search should finish');
               flag=1;
           end
        end
    end

if flag ==1
   break 
end  

    end

if flag ==1
   break 
end
end


 U1=otimes(maxpower(A,m),otimes(W,maxpower(B,n)));
U;
U1-U;
 