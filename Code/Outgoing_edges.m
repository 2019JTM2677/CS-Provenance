function [array] = Outgoing_edges(node_b,n)
    

    if node_b==n
        fprintf("No outgoing edges from destination\n")
        return
    end
    array=[];
    for k=1:n-1
       if k<node_b
           array = [array (k-1)*(n-1)+ node_b-1];
       elseif k==node_b
           array = [array (n-1)*node_b];
       else 
           array = [array (k-1)*(n-1)+ node_b];
       end
    end

end