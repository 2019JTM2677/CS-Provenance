% Function for OMP algorithm 
% K sparsity, x recovered signal, y measured signal, A matrix
function [x] = Path_aware_OMP (K,y,A,n,alpha)

[~,E] = size (A) ;
x = zeros (E,1) ;
Residue = [];
B=[];
kk = []; % Support set

% Iterating K times 
for J = 1 : K
    %Index Search
    if J==1
        [~ ,index] = max(abs(A' * y)) ;
    else
        %[~,index] = max(abs(A' * Residue(:,J-1)));
       
        if J==2
            g11 = []; g12=[];g31=[];g32=[];
            edges=[];reverse_group=[];
        end
        [index,g11,g12,g31,g32,reverse_group,edges] = max_index(A,Residue(:,J-1),g11,g12,g31,g32,reverse_group,edges,kk(J-1),n,alpha);
        %fprintf("Final: ");display(edges);
        %}
    end
        
        kk (J) = index;
        %kk = union(index,kk);
    
    %Residue Update
    B = A(:,kk);
   
    x_cap = B\y;%pinv(B)*y;      % expected x

    Residue(:,J) = y - B*x_cap;  % residue 
    
end
% Final x after recovery
x(kk)=x_cap;
end