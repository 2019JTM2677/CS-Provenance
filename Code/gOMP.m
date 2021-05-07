% Generalized OMP
% https://arxiv.org/pdf/1111.6664.pdf

function [x]=gOMP(K,y,A,N)
%A: Sensing matrix
%y: measurement vector
%K: Sparsity
%N: No. of indices chosen in each iteration
[m,E] = size (A) ;
x = zeros (E,1) ;
if (N > K) || (N > m/K)
 fprintf("Give correct value for N");
 return 
end

k=0;
Residue =y;
B=[];      
lambda=[];% support set
e = 10^-5;

while k< min(K,m/N) && norm(Residue) > e
    k=k+1;
    [~ ,index] = maxk(abs(A' * Residue),N);
    
    lambda = union(lambda,index);
    
    B= A(:,lambda);
    %display(B);
    %B_mod = B'*B;
    %y_mod = B'*y;
    %x_cap = B_mod\y_mod;
    x_cap = B\y;%pinv(B)*y;
    Residue = y - B*x_cap;
    
end
x(lambda)=x_cap;
end