function A=GenerateRandomPolynomial(n,range1,range2)
% GenerateRandomMatrix generates an mxn matrix with entries between
%range 1 and range 2
A=randi([range1, range2],1,n);