function [array] = Incoming_edges(node_a,n)

    array=[];
    if node_a ==n
        array = (n-1):(n-1):(n-1)^2;
    else 
        array = 1+ (node_a - 1)*(n-1): node_a*(n-1) -1;
    end    

end