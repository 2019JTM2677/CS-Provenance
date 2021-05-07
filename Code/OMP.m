% Function for OMP algorithm 
% K sparsity, x recovered signal, y measured signal, A matrix
function [x] = OMP (K,y,A)

[m,n] = size (A) ;
x = zeros (n,1) ;
Residue = [];%zeros(m,K);
B=[];
%x_cap=zeros(K,1);
kk = [];

% Iterating K times 
for J = 1 : K
    
    %Index Search
    if J==1
        [~ ,index] = max(abs(A' * y)) ;
    else
        [~,index] = max(abs(A' * Residue(:,J-1)));
    end
        %kk (J) = index ;
        kk = union(index,kk);
    
    %Residue Update
    %B = [B A(:,kk(J))];          % Basis matrix
    B = A(:,kk);
    
    B_mod = B'*B;
    y_mod = B'*y;
    % x_cap = B_mod\y_mod;%B'*y  ;         % x^ in each iteration
    x_cap = B\y;%pinv(B)*y;      % expected x

    Residue(:,J) = y - B*x_cap;%new_x';  % residue 
    
end
% Final x after recovery
x(kk)=x_cap;
end