% PL-OMP for double edge embedding
function [x_de] = PL_OMP_DE(n,h_de,h,y,A,L,B_mat,C_mat,A_1,A_2,temp_a,temp_b,temp_c,temp_d,temp_e,src_id,dest_id,penultimate_id)

    x = OMP_modified(h_de-1,y,A,L);
    x(abs(x)<=0.001)=0;
    x(abs(x)>0.001)=1;
    
    jj=[];    % Invalid path indices from set of solutions
    x_de=[]; count=0;
    
    % Missing double edge contraints
    for j=1:size(x,2)
        x_n_sq = DE_to_SE(x(:,j),n,h,penultimate_id,temp_a,temp_b,temp_c,temp_d,temp_e,A_1,A_2); 
        [opt,src_connect_node, dest_connect_node,nodes_in_path] = one_missing_link(x_n_sq,n,h,src_id,dest_id,B_mat,C_mat,'DE'); 
        
        nodes = 1:n;
        nodes_left = setdiff(nodes, nodes_in_path);      % options for middle node
        if opt==1
            missing_index = zeros(size(nodes_left,2),1); % options for missing double edge index
            for ni = 1:size(nodes_left,2)
                mid_node = nodes_left(ni);
                if src_connect_node < mid_node
                    strt = (src_connect_node - 1)*(n-2);
                elseif src_connect_node > mid_node
                    strt = (src_connect_node - 2)*(n-2);
                end
                
                if dest_connect_node < mid_node && dest_connect_node < src_connect_node
                    temp_index = dest_connect_node;
                elseif dest_connect_node > mid_node && dest_connect_node > src_connect_node
                    temp_index = dest_connect_node - 2;
                else
                    temp_index = dest_connect_node -1;
                end
                missing_index(ni) = (mid_node-1)*(n-2)^2 + strt + temp_index;
                count = count +1;
                x_de(:,count) = x(:,j);
                x_de(missing_index(ni),count)=1;
            end
            
        else
            jj=[jj j];
        end
    end
    
    x(:,jj)=[];

    %Arrange in ascending residual norm
    residual_norm=zeros(1,size(x_de,2));
    for j=1:size(x_de,2)
        final_residue = y - A*x_de(:,j);
        residual_norm(j) = norm(final_residue);
    end
    [~,order]= sort(residual_norm);
    x_final = x_de;
    x_de = x_final(:,order);
end

                
