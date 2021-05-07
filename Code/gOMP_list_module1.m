function [x,ratio_bad_by_total] = gOMP_list_module1(K,y,A,N,n,temp_e,temp_f,temp_g,temp_dest,B_mat,C_mat,src_id)
% PL-gOMP
% least square & residue calculation at last iteration
%A: Sensing matrix
%y: measurement vector
%K: Sparsity
%N: No. of indices chosen in each iteration
[m,E] = size (A) ;
%x = zeros (n,1) ;
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

while k < tot_no_of_iter
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
                
                index_combo=nchoosek(new_index(:,i),N);
                for j=1:L
                    indices{end+1}=union(lambda{i},index_combo(j,:));    
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
    
    % Path constraints check at last iteration
    if k == tot_no_of_iter
        jj=[];final_lambda={};
        size_before = size(lambda,2);
        x = zeros (E,size(lambda,2)) ;
        %fprintf("before:");display(size(lambda,2));
        for i=1:size(lambda,2)
           %display(lambda{i})
           x(lambda{i},i)= 1;
           x_nminus1_sq = x(:,i);
           x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
           %Adjacency_constraint(x_n_sq,B_mat,C_mat,K,src_id,'gomp_module_1');
           A_mat = diag(x_n_sq);
           Z = B_mat*A_mat*C_mat;
           Za = Z^K;
           if Za(end,src_id)==0
                jj=[jj i];
           else
               final_lambda{end+1}=lambda{i};
           end
        end
        x(:,jj)=[];
        
        lambda = final_lambda;
    end
    B=[];
    x_cap=[];x_cap_save={};
    for i=1:size(lambda,2)
        %display(k)
        %display(lambda{i})
        B= A(:,lambda{i});
        x_cap = B\y;%; pinv(B)*y   %Estimation
        x_cap_save{i}=x_cap;
        Residue(:,i) = y - B*x_cap;     %Residual update
        %display(Residue)
    end
end
size_after = size(lambda,2);
ratio_good_by_total = size_after/size_before;
ratio_bad_by_total = 1 -ratio_good_by_total; % fraction of bad solutions pruned in path constraints list 
x = zeros (E,size(lambda,2)) ;

for i=1:size(lambda,2)
    x(lambda{i},i)= x_cap_save{i};
end

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