% Stage wise orthogonal matching pursuit
% http://sparselab.stanford.edu/SparseLab_files/local_files/StOMP.pdf

function [x]= StOMP(K,y,A)


    [m,n] = size (A) ;
    x = zeros (n,1) ;
    x_cap=[];
    Residue = y;
    curr_I=[];   % support set
    S = 10; % fixed No. of iterations, if S=m, then StOMP OMP gives same x
    s=0;e=10^(-5);
    rho =K/m;
    delta=m/n;
    alpha_0 = delta*(1-rho)/S;
    q = min((m-K)/K,0.5);
    
    for i=1:S
       %s=s+1;
       val =  sqrt(m) .* abs(A' * Residue) ./norm(Residue)  ;      % projection values
       
       thr = norminv(1 - alpha_0/2, 0, 1);                  %Threshold parameter
       %sigma = norm(Residue)/sqrt(m);          % Formal noise level
       %thr = t;%*sigma;
       %thr= fdrthresh(val, q);
       % Apply hard thresholding
        I = find(val > thr);
    
      
       J = union(curr_I,I);       % Hard thresholding
       
       if (length(J) == length(curr_I)) 
            %done 
            break
       else
            curr_I = J;
       
       %}
       display(I);
       %x(I) = pinv(A(:,I))*y;
       %Residue = y-A*x;
       %x_cap = pinv(A(:,curr_I))*y;
       x_cap = A(:,curr_I)\y;
       Residue = y - A(:,curr_I)*x_cap;
       end
       if norm(Residue)<=e
           break
       end
    end
x(curr_I)=x_cap;
end
