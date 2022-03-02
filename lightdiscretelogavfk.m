function power=lightdiscretelogavfk(A,V,F)
% Finds k such that A=V F^k
% works when lambda(F) is not equal to 0
% currently only for finite entries

%eps=10^{-9}
eps=0.000000001;

%disp('Discretelog starts here.');
%A
%V
%F
[m,n]=size(F);
if m~=n
    error('The third input argument should be square matrix');
end

[m,n1]=size(V);
if n1~=n
    error('The dimensions of the second and third input argument dont match');
end

[m2,n2]=size(A);
if m2~=m || n2~=n
 error('The dimensions of the input arguments do not match');
end


disp('applying policy iteration');
[chi,x,criticalcycle]=policyIteration(F)
lambda=max(chi);
disp('policy iteration finished');
if abs(lambda)<eps
    error('The mcm of the third argument is equal to 0. The method does not work in this case.');
end

I=maxeins(n);
FF=F-lambda;
% Finding an eigenvector of FF
% disp('Now finding an eigenvector of FF');
% M=metricmatrix(FF);
% 
% for i=1:n
% if abs(M(i,i))<eps;
%     break
% end
% end
% 
% FFstar=max(I,M);
% eigvector=FFstar(:,i);


%disp('Now finding a critical cycle in FF');
%cycle =critical(FF,eigvector);
%cycle =critcircuit(FF)
[~,l]=size(criticalcycle)
c=l+1
criticalcycle(l+1)=criticalcycle(1)
zz=zeros(n,1)
for i=1:c
    zz(criticalcycle(i))=1
end
%zz

disp('Now using the CSR to catch the power');
disp('Computing C and S');

%U=maxfloyd(powmaxplus(FF,c));
U=maxfloyd(maxpower(FF,c-1))
C=column(U,zz);


S=-inf*ones(n,n);
 for i=1:c-1
S(criticalcycle(i),criticalcycle(i+1))=FF(criticalcycle(i),criticalcycle(i+1));
 end
%S
for j=1:c-1
Y(:,j)=A(:,criticalcycle(j));
end
% Now we will compare a column with index in the 
% critical cycle ("cycle") of Areduced with V otimes C otimes S^d
% for all possible values of d. 
% When they are different by the scalar, we deduce the power from the 
% difference.
disp('Now catching the power');
c
newflag=0;
 for d=1:c-1   
 W=otimes(C,binmaxpower(S,d))
 Acand=otimes(V,W)
 %Y
 for j=1:c-1
 X(:,j)=Acand(:,criticalcycle(j))
 end
 %X
 
 Z=Y-X
 % will not be good if we have -inf entries..
 Mu1=Z(1,1)
 d
 Z=Z-Mu1
 if norm(Z)<eps
     newflag=1;
     break
 end
 end


if newflag==0
    power=-1;
else
power=round(Mu1/lambda);
end

end



 