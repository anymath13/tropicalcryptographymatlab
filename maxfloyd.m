function g=maxfloyd(A)
% Computing the max-plus star of a given matrix or 
% detecting a positive cycle.

eps=0.0000000001;
[n,n]=size(A);

I=zeros(n,n);

for i=1:n
    for j=1:n
        if i==j
             I(i,j)=0;
        else
             I(i,j)=-inf;
        end
    end
end

g=A;
for k=1:n
    for i=1:n
        for j=1:n
            if (i~=k)&&(j~=k)
                if g(i,k)+g(k,j)>g(i,j)
                    g(i,j)=g(i,k)+g(k,j);
                end
            end
            if i==j && g(i,i)>eps
            error('Terminated after detecting a positive cycle');
            end
        end
    end
end
g=max(I,g);

%G=B;
%g=B;
%for i=1:n
%    G=otimes(B,G);
%    g=max(G,g);
%end
%g=max(I,g);
%kleen=max(I,gamma);
%disp('metricmatrix');
%disp(gamma);



