function [path_arr] = path_array(path,E,n)
    
    path_arr = zeros(E,1);
    Ap = zeros(n-1,n-1);
    for i=1:n-1
        for j=2:length(path)
            curr_node = path(j);
            prev_node = path(j-1);
            if curr_node == n
                Ap(prev_node,n-1) = 1;
            else
                if prev_node < curr_node
                    index = prev_node;
                else
                    index = prev_node-1;
                end
                Ap(curr_node,index) = 1;
            end
        end
    end
    path_arr=reshape(Ap',[E,1]);
end