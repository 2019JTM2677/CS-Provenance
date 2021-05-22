% Function for PL-OMP with h-2 iterations
function [x_final] = PL_OMP2(n,h,y,A,L,B_mat,C_mat,temp_e,temp_f,temp_g,temp_dest,src_id,dest_id,method)
 
    x = OMP_modified(h-2,y,A,L);
    x(abs(x)<=0.001)=0;
    x(abs(x)>0.001)=1;
    x_final=[];count=0;
    
    jj=[];
    for j=1:size(x,2)
        x_nminus1_sq = x(:,j);
        x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
        [opt,src_connect_node, dest_connect_node,n1,n2,nodes_in_path] = two_missing_link(x_n_sq,n,h,src_id,dest_id,B_mat,C_mat,method);
        if opt==1
            nodes = 1:n;
            nodes_left = setdiff(nodes, nodes_in_path);      % options for middle node
            % Find missing index src_conn_node->mid & mid->dest
            for ni = 1:size(nodes_left,2)
                mid_node = nodes_left(ni);
                first_missing_index = FindIndex(n,src_connect_node,mid_node,dest_id);
                second_missing_index = FindIndex(n,mid_node,dest_connect_node,dest_id);
                count = count +1;
                x_final(:,count) = x(:,j);
                x_final([first_missing_index,second_missing_index],count)=1;
            end
            
        elseif opt==2
            first_missing_index = FindIndex(n,src_connect_node,n1,dest_id);
            second_missing_index = FindIndex(n,n2,dest_connect_node,dest_id);
            if first_missing_index > (n-1)^2 || second_missing_index > (n-1)^2
                fprintf("Error in index src_conn:%d,n1:%d, dest_con:%d,n2:%d\n",src_connect_node,n1,dest_connect_node,n2)
                fprintf("Links in x:");display(find(x(:,j)))
            end
            count = count +1;
            try 
                x_final(:,count) = x(:,j);
            catch
                fprintf("ERror in x\n")
                display(size(x(:,j)));
                display(size(x_final));
            end
            x_final([first_missing_index,second_missing_index],count)=1;
        else
            jj=[jj j];
        end
        
    end
    
    x(:,jj)=[];

    %Arrange in ascending residual norm
    residual_norm=zeros(1,size(x_final,2));
    for j=1:size(x_final,2)
        
        try
            final_residue = y - A*x_final(:,j);
        catch
            fprintf("ERRORRR\n")
            display(size(x_final(:,j)))
            display(size(A))
        end
        residual_norm(j) = norm(final_residue);
    end
    [~,order]= sort(residual_norm);
    x_final = x_final(:,order);
end