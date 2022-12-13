function C=otimes(A,B)

[m,n]=size(A);
[k,l]=size(B);
     Ax = [];     
    Ax1 = []; 
for i=1:m
    for j=1:l
        Ax1 = A(i,:) + B(:,j)';
        Ax(i,j) = max(Ax1);
        
    end
end  
C=Ax;

