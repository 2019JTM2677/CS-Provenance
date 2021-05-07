function [x]=gOMP_list_module2(K,y,A,N)
% L-gOMP
%A: Sensing matrix
%y: measurement vector
%K: Sparsity
%N: No. of indices chosen in each iteration
[m,E] = size (A) ;
if (N > K) || (N > m/K)
 fprintf("Give correct value for N");
 return 
end

k=0;
Residue =[];
B=[];      
lambda={};% support set
e = 10^-5;
N_d=N+1;
L=nchoosek(N_d,N);
tot_no_of_iter = min(K,m/N);

while k< tot_no_of_iter
    k=k+1;
    new_index=[];index_combo=[]; indices={};
    if k==1
        [~ ,new_index] = maxk(abs(A' * y),N_d);
        index_combo = nchoosek(new_index,N);
        for j=1:L
            indices{end+1}=index_combo(j,:);
            %display(indices{end})
        end
    else
        for i=1:size(lambda,2)
            if norm(Residue(:,i)) > e
                [~ ,new_index(:,i)] = maxk(abs(A' * Residue(:,i)),N_d);
                %display(abs(A' * Residue(:,i)));
                %display(val)
                %display(index(:,i))
                %index_combo=[repmat(lambda{i},L,1) nchoosek(new_index(:,i),N)]; %Augmentation
                index_combo=nchoosek(new_index(:,i),N);
                for j=1:L
                    %indices{end+1}=index_combo(j,:);
                    indices{end+1}=union(lambda{i},index_combo(j,:));
                    %display(indices{end})
                end
            else 
                index_combo = lambda{i};
                indices{end+1}=index_combo;
                %display(indices{end})
            end
            %display(d); display(size(indices));
        end
    end
    %fprintf("iter: %d\n",k);
    lambda = indices;
    
    B=[];
    x_cap=[];x_cap_save={};
    
    for i=1:size(lambda,2)
        B= A(:,lambda{i});
        x_cap = B\y;%; pinv(B)*y   %Estimation
        x_cap_save{i}=x_cap;
        Residue(:,i) = y - B*x_cap;     %Residual update
    end
end

x = zeros (E,size(lambda,2)) ;
for i=1:size(lambda,2)
    x(lambda{i},i)= x_cap_save{i};
end

%display(residual_norm)
%{
x(abs(x)<=0.001)=0;
x(abs(x)>0.001)=1;
%display(residual_norm)
for j=1:size(x,2)
    final_residue(:,j) = y - A*x(:,j);
    residual_norm(j) = norm(final_residue(:,j));
end

[~,b]= sort(residual_norm);
x_final = x;
x = x_final(:,b);
%}
end