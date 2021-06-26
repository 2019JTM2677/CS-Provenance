
function index = reverse_edge(node_a,node_b,n)

    if node_b < node_a
        i_const = node_b;
    elseif node_b > node_a
        i_const = node_b -1;
    end
    % Index for node_b->node_a
    %fprintf("%d %d\n",node_a,node_b)
    try
        index = (n-1)*(node_a -1)+i_const;
    catch
        fprintf("%d %d\n",node_a,node_b)
    end
end