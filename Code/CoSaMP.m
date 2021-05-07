% CoSaMPs
% http://www.cmc.edu/pages/faculty/DNeedell
function [x] = CoSaMP(K,y,A)

    [m,n] = size (A) ;
    x = zeros (n,1) ;
    Residue = y;
    I=[];   % support set
    %support=[];
    kk=[];
    e = 10^-5;
    k=0;
    %for i=1:4
    while (norm(Residue) > e) && (k<10)
        k=k+1;
       [~ ,index] = maxk(abs(A' * Residue),2*K); 
       display(index);
       I = union(index,kk);
       display(I);
       x_cap = A(:,I)\y;%pinv(A(:,I))*y;
       display(x_cap);
       [~,support] = maxk(abs(x_cap),K);
       display(support);
       kk = I(support);
       x_cap_pruned = x_cap(support);
       display(kk)
       %display(x_cap_pruned)
       Residue = y - A(:,kk)*x_cap_pruned;
       %fprintf("res: %d\n",norm(Residue))
    end
x(kk)=x_cap_pruned;
end