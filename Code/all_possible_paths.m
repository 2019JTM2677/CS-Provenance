function path_matrix = all_possible_paths(n,h,src_id,dest_id)
    path_matrix=[];
    node_id = 1:n;
    node_combo = nchoosek(node_id,h+1);
    %display(size(node_combo,1))
    for i=1:size(node_combo,1)
        path_matrix = [path_matrix;perms(node_combo(i,:))];
    end
    
    path_matrix(path_matrix(:,1)~=src_id,:)=[];
    path_matrix(path_matrix(:,end)~=dest_id,:)=[];

end