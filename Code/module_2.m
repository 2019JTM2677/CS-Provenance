% Function for L-OMP
function [x,size_v2 ] = module_2(n,h,y,Ar,L,B,C,temp_e,temp_f,temp_g,temp_dest,src_id)

    x_OMP_mod = OMP_modified(h,y,Ar,L);
    x_OMP_mod(abs(x_OMP_mod)<=0.001)=0;
    x_OMP_mod(abs(x_OMP_mod)>0.001)=1;
    %fprintf("List of solutions:")
    %display(x_OMP_mod)
    size_v2 = size(x_OMP_mod,2);
    
    jj=[];        % Invalid path indices 
    
    % Path constraints 
    for j=1:size(x_OMP_mod,2)
        x_nminus1_sq = x_OMP_mod(:,j);
        x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);

        opt = Adjacency_constraint(x_n_sq,B,C,h,src_id);
        if opt==1
           %fprintf("Valid h len path");
        else
            %fprintf("Invalid");
            jj=[jj j];
        end
    end
    x_OMP_mod(:,jj)=[];         % Remove invalid solutions
    
    % Arrange solutions in ascending residual norm 
    residual_norm=zeros(1,size(x_OMP_mod,2));
    for j=1:size(x_OMP_mod,2)
        final_residue = y - Ar*x_OMP_mod(:,j);  %Use this only for solutions that verify path constraints
        residual_norm(j) = norm(final_residue);
    end
   
    [~,order]= sort(residual_norm);
    x = x_OMP_mod(:,order);
    
end