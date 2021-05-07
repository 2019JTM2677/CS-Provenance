classdef node
   properties
      Node_id {}
      Edge_id {}
      reachable {}
      curr_neighbour {}
     
   end
   methods
       function obj = node(val1,val2,val3)
            if nargin == 3
                obj.Node_id = val1;
                obj.Edge_id = val2;
                obj.reachable = val3;
                %obj.curr_neighbour = val4;
            end
            obj.curr_neighbour = randi(1,sum(obj.reachable));
        end
       function r = find_reachable(obj)
        for i = 1: length(obj.reachable)
            if obj.reachable(i)==1
                r = [r i];
            end
        end
      end
      
   end
end