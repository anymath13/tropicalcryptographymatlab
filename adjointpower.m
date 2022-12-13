function A = adjointpower(H,k)
%kth adjoint power of a matrix (H oplus H^2 oplus ... oplus H^k)
[m,n]=size(H);
if m~=n
    error('Matrix should be square');
end

if k<1
    error('Only positive powers can be taken')
end

A=H;

for i=2:k
Apower=otimes(A,H);
A=max(A,Apower);
end


    