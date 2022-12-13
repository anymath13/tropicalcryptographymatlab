function C=add(A,B)

[m,n]=size(A);
Ax=[];
for i=1:m
    for j=1:n
        Ax(i,j)=max(A(i,j),B(i,j));
    end
end
C=Ax;