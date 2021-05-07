function [opt,src_conn_id,dest_conn_id,nodes_in_path] = one_missing_link(x_n_sq,n,h,src_id,dest_id,B,C,method)
    
    opt=0;%r_id=[];c_id=[];
    dest_conn_id=[];
    src_conn_id=[];
    nodes_in_path=[];    %EE = size(x_n_sq);
    
    A = diag(x_n_sq);
    Z = B*A*C;
    %display(Z);
    
    if any(sum(Z)>1) || any(sum(Z,2)>1)        % More than one non zero in any column of Z,i.e. more than 1 outgoing edges from same node
        %fprintf("invalid path no need for further check")
        %fprintf("BUG")
        
    else
    
        dest_mat=[]; 
        src_mat=[];
        switch lower(method)
            case 'se'
                param=h-1;
            case 'de'
                param=h-2;
        end
        for nd=1:param
            Za = Z^nd;
            dest_mat=[dest_mat;Za(end,:)];
            src_mat=[src_mat Za(:,src_id)];
        end
        %display(dest_mat);
        %display(src_mat);
        
        r = zeros(1,param);
        r(1)=1;
        %display(r);
        c = zeros(param,1);
        c(1)=1;
        rk=0;ck=0;
        %r_id=[];c_id=[];
        %nodes_in_path=[]; % all nodes appeared in path
        nodes_in_path = union(nodes_in_path,[src_id,dest_id]);
        src_flag=0;
        dest_flag=0;
        for k=1:param
            if k==1 || src_flag==1
                I=ismember(src_mat,r,'rows');
                %display(I)
                if  sum(I)==1 
                    [~,r_id]=ismember(r,src_mat,'rows'); % node connected to src
                    nodes_in_path = union(nodes_in_path,r_id);
                    if r_id~=src_id && r_id ~= dest_id  % node connected to src in h-1 hops can't be src or dest
                        src_conn_id = r_id;
                        src_flag=1;
                        rk=rk+1;
                        r=circshift(r,1);
                    else
                        src_flag=0;
                    end
                else
                    src_flag=0;
                % If no match found in first iteration then we don't need to go for next iteration as if hop len 1 doesnt exist 2,3 can't exist
                % similarily if 2 doesn't exist 3  and onwards can't
                end
            end

            if k==1 || dest_flag==1
            J=ismember(dest_mat',c','rows') ;
                if sum(J)==1
                    [~,c_id]=ismember(c',dest_mat','rows'); % node connected to dest
                    nodes_in_path = union(nodes_in_path,c_id);
                    if c_id~=src_id && c_id ~= dest_id
                        dest_conn_id = c_id;
                        dest_flag=1;
                        ck=ck+1;
                        c=circshift(c,1);
                    else
                        dest_flag=0;
                    end
                else
                    dest_flag=0;
                end
            end


        end
        if isempty(src_conn_id)
            src_conn_id = src_id;
        end
        if isempty(dest_conn_id)
            dest_conn_id = dest_id;
        end


        switch lower(method)
        case 'se'
            if rk+ck==param && src_conn_id~=dest_conn_id
                %display(r)
                %display(c)
                opt=1;
            end
        case 'de'
            if rk+ck==param
                %display(rk)
                %display(ck)
                opt=1;
            end
        end
    end
end