function [packet]= edge_embed(m,path,packet,nodes) 
    packet.provenance = double(zeros(m,1)); % m= no. of row in A
    for i=2:length(path)
        if i==length(path)
            packet.provenance = packet.provenance + double(nodes(path(i)).Edge_id(:,(length(nodes)-1)*nodes(path(i-1)).Node_id)) ;
        else
            update_prov = calculate_index(nodes(path(i-1)), nodes(path(i)));
            %disp(update_prov)
            %disp(packet.provenance)
            packet.provenance = packet.provenance +  double(update_prov);
        end
        
       
    end
   % fprintf("Final provenance\n");disp(packet.provenance)
end

function [edge] = calculate_index(prev_node,curr_node)
    Edge = curr_node.Edge_id;
    if prev_node.Node_id < curr_node.Node_id
        index = prev_node.Node_id;
    else
        index = prev_node.Node_id -1;
    end
    edge = Edge(:,index);
end