function [opt,src_conn_id,dest_conn_id,n1,n2,nodes_in_path] = two_missing_link(x_n_sq,n,h,src_id,dest_id,B,C,method)
    
    opt=0; % If opt=1, brute force mid node. If opt=2, connect src_conn_id->n1 & n2->dest_conn_id
    dest_conn_id=[];
    src_conn_id=[];
    nodes_in_path=[];    
    n1=0;n2=0;
    
    A = diag(x_n_sq);
    Z = B*A*C;
    %display(Z);
    
    if any(sum(Z)>1) || any(sum(Z,2)>1) || nnz(Z)~=h-2  % More than one non zero in any column of Z,i.e. more than 1 outgoing edges from same node
        %fprintf("invalid path no need for further check")
    else
        dest_mat=[]; 
        src_mat=[];
        
        % Set parameter for missing-link
        switch lower(method)
            case 'se'
                param=h-2;
        end
        Zb=[];
        for nd=1:param
            Za = Z^nd;
            Zb = cat(3,Zb,Za);
            dest_mat=[dest_mat;Za(end,:)];
            src_mat=[src_mat Za(:,src_id)];
        end
        %display(dest_mat);
        %display(src_mat);
        
        r = zeros(1,param);
        r(1)=1;
        c = zeros(param,1);
        c(1)=1;
        rk=0;ck=0;
        
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
                    if r_id~=src_id && r_id ~= dest_id   % node connected to src in h-1 hops can't be src or dest
                        src_conn_id = r_id;
                        src_flag=1;
                        rk=rk+1;
                        r=circshift(r,1);
                    else
                        src_flag=0;
                    end
                else
                    src_flag=0;
               
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
                opt=1;
                %fprintf("Consecutive missing edges\n")
            elseif rk+ck<param && src_conn_id~=dest_conn_id
                %display(find(sum(Z)));
                %display(nodes_in_path)
                %fprintf("Separate missing edges\n")
                nodes_left = setdiff(find(sum(Z)),nodes_in_path);
                %display(nodes_left);
                P=Zb(:,nodes_left,:);
                d=zeros(h-2,1);
                d(h-2-(rk+ck))=1;
                pattern(1,1,:)=d;
                whereToFind=all(P==pattern,3);
                if nnz(whereToFind)==1
                    %fprintf("Found\n")
                    [n2 n1]=find(whereToFind);
                    n1 = nodes_left(n1);
                    if length(unique([n1,n2,src_conn_id,dest_conn_id]))==4
                        opt=2;
                    end
                end
                %display(n1);
                %display(n2);
                %fprintf("nodes: %d %d",n1,n2);
                %fprintf("Size of paths: %d %d %d",rk,h-2-rk-ck,ck);
            end
        
        end
    end
end