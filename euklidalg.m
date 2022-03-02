function [d,p,q,r,s] = euklidalg(a,b,c);
% Solution of Diophantine equation ax + by = c
% by means of Euclidean algorithm
% The greatest common divisor is the last nonzero remainder
% Back substitution:
% coefficients of the linear combination p, q
% r = -b/d, s = a/d
%
% test of equation correctness:
if isempty(find(a)) | isempty(find(b))
display('Diophantine equation is not given correctly!')
p=0; q=0; d=0; s=0; r=0;
return
end
alfa=a; beta=b; rem=a;
% search for the greatest common divisor:
i=0;
n=abs(length(a)-length(b))+1;
if n == 1
n=n+1;
end
% test of zero remainder:
while norm(rem,inf) > 100*eps
i=i+1;
% elimination the zero leading coefficients:
ind=find(beta);
beta=beta(ind(1):length(beta));
% quotient and remainder:
[quot,rem]=deconv(alfa,beta);
i0=1+n-length(quot);
% storing of quotients:
qq(i,i0:n)=quot;
% shift of polynomials:
alfa=beta; beta=rem;
end
% recurrent computation of the coefficients
d=alfa; p=0; q=1; m=i-1
for i=m:-1:1
r=p; p=q
% formal rearrangement for polynomial sum executing:
rr=zeros(1,length(qq(i,:))+length(p)-1);
rr(length(rr)-length(r)+1:length(rr))=r;
% computation of further element of the sequence:
q=rr-conv(qq(i,:),p);
end
% normalization of polynomial:
ind=find(q); q=q(ind(1):length(q));
% general solution of the reduced equation:
r=-deconv(b,d); s=deconv(a,d);
return