% Function for double edge embedding
function [packet]= double_edge_embed(m,path,packet,nodes) 
    packet.provenance = double(zeros(m,1)); % m= no. of row in A
    flag = 1;
    n = size(nodes,2);
    
    for i=2:length(path)-1
        %if i==length(path)
            %packet.provenance = packet.provenance + double(nodes(path(i)).Edge_id(:,(length(nodes)-1)*nodes(path(i-1)).Node_id)) ;
        %else
            %fprintf("curr node:%d",path(i));
            %fprintf(" prev node:%d",path(i-1));
            %update_prov = calculate_index(nodes(path(i-1)), nodes(path(i)),nodes(path(i+1)));
            
            %disp(update_prov)
            %disp(packet.provenance)
            
            if flag
                update_prov = calculate_index(nodes(path(i-1)), nodes(path(i)),nodes(path(i+1)),n);
                packet.provenance = packet.provenance +  double(update_prov);
            end
            
            flag = ~flag;
            if rem(length(path),2)==0  % odd hop length
               if path(i+1)==n
                   flag=0;
               end
            end
        %end
        
       
    end
   % fprintf("Final provenance\n");disp(packet.provenance)
end

function [edge] = calculate_index(prev_node,curr_node, next_node,n)
    Edge = curr_node.Edge_id;
    
    if prev_node.Node_id < curr_node.Node_id 
        strt = 1+(prev_node.Node_id - 1)*(n-2);
        %finish = start + (n-3); 
    elseif prev_node.Node_id > curr_node.Node_id
        strt = 1+(prev_node.Node_id - 2)*(n-2);
        %finish = start + (n-3);
        %index = prev_node.Node_id -1;
        %fprintf("Error in path");
    end
    
    if next_node.Node_id < curr_node.Node_id && next_node.Node_id < prev_node.Node_id
        index = next_node.Node_id;
    elseif next_node.Node_id > curr_node.Node_id && next_node.Node_id > prev_node.Node_id
        index = next_node.Node_id - 2;
    else
        index = next_node.Node_id -1;
    end
    
    final_index = strt + index -1;
    edge = Edge(:,final_index);
    
end