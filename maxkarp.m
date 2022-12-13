function lambda=maxkarp(A)
% Syntax lambda=maxkarp(A)
% Computes eigenvalue in min-plus algebra using Karp's algorithm
[n,m]=size(A);
if (n~=m)
    error('Matrix must be square')
end
B=A
t=-inf*ones(n,1)
for i=1:n
    C(i,:,:)=B;
    B=otimes(A,B)
end

t(1)=C(n,1,1)/n
for i=1:n
    for k=1:n-1
        t(i)=min(t(i), (C(n,i,1)-C(k,i,1))/(n-k));
    end
end
lambda=max(t);