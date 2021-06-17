% Function for modified OMP algorithm 
% h sparsity, x recovered signal, y measured signal, A matrix
% L no. of max projection indices chosen for each projection
function [x] = OMP_modified (h,y,A,L)

    [~,E] = size (A) ;
    x = zeros (E,L^h) ;

    Residue = [];
    S = [];        % Support set for list of solutions

    for j = 1 : h
        x_cap=[];
        ind=[];
        %fprintf("iter:%d",j)
        % Projection 
        if j==1
            [~ ,index] = maxk(abs(A' * y),L); %choose first L max elements
            %[c,d]=sort(abs(A' * y))
        else
            for l=1:L^(j-1)
                [~,ind(:,l)] = maxk(abs(A' * Residue(:,l,j-1)),L);
                %[c,d]=sort(abs(A' * Residue(:,l,j-1)))
            end
            ind = reshape(ind,[L^j,1]);
            S = repelem(S,L,1);
            index = ind;
        end
        S (:,j) = index  ;  % Add the newly selected indices to support set
   
        %display(S);
        
        for l = 1:L^j
            %C = A(:,kk(l,:))'*A(:,kk(l,:));%C(:,:,l)'*C(:,:,l);
            %y_mod = A(:,kk(l,:))'*y;
            %x_cap(:,l) = C\y_mod;%B'*y  ;         
            x_cap(:,l) = A(:,S(l,:))\y; %pinv(C(:,:,l))*y;
            Residue(:,l,j) = y - A(:,S(l,:))*(x_cap(:,l));%new_x';  % residue 
        end
    end
    
    for l=1:L^h
        for i=1:size(S,2)
            x(S(l,i),l)= x_cap(i,l);     % Final x after recovery
        end
    end
    %{
    final_residue=zeros(m,size(x,2));
    residual_norm=zeros(1,L^h);
    for j=1:size(x,2)
        final_residue(:,j) = y - A*x(:,j);  %Use this only for solutions that verify path constraints
        residual_norm(j) = norm(final_residue(:,j));
    end
    %display(residual_norm)
    [~,b]= sort(residual_norm);
    x_final = x;
    x = x_final(:,b);
    %}
end