% Function for PL-OMP
function [x] = module_1(n,h,y,A,L,B_mat,C_mat,temp_e,temp_f,temp_g,temp_dest,src_id,dest_id,method)
 
    x = OMP_modified(h-1,y,A,L);
    x(abs(x)<=0.001)=0;
    x(abs(x)>0.001)=1;
    
    jj=[];
    for j=1:size(x,2)
        x_nminus1_sq = x(:,j);
        x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
        [opt,src_connect_node, dest_connect_node, ~] = one_missing_link(x_n_sq,n,h,src_id,dest_id,B_mat,C_mat,method);
        if opt==0
            jj=[jj j];
            
        else
            if dest_connect_node == dest_id
                missing_index = (n-1)*src_connect_node;
            elseif src_connect_node < dest_connect_node
                missing_index = (n-1)*(dest_connect_node-1)+src_connect_node;
            else 
                missing_index = (n-1)*(dest_connect_node-1)+src_connect_node-1;
            end
            try
                x(missing_index,j)=1;
            catch
                fprintf("mod 1:")
                display(missing_index)
            end
        end
        
    end
    
    x(:,jj)=[];
    residual_norm=zeros(1,size(x,2));
    for j=1:size(x,2)
        final_residue = y - A*x(:,j);
        residual_norm(j) = norm(final_residue);
    end
    [~,order]= sort(residual_norm);
    x_final = x;
    x = x_final(:,order);

end