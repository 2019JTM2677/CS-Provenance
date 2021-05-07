function [path_arr] = path_array_de(path,E,n)
    path_arr = zeros(E,1);Ap = zeros((n-2)^2,n-1);
    flag = 1;
    %for i=1:n-1
        
        for j=2:length(path)-1
            
            curr_node = path(j);
            prev_node = path(j-1);
            next_node = path(j+1);
            
            if flag
                %fprintf("Embedding nodes")
                %disp(curr_node);
                if prev_node < curr_node
                    strt = 1+(prev_node - 1)*(n-2);
                elseif prev_node > curr_node
                    strt = 1+(prev_node - 2)*(n-2);
                end
                
                if next_node < curr_node && next_node < prev_node
                    index = next_node;
                elseif next_node > curr_node && next_node > prev_node
                    index = next_node - 2;
                else
                    index = next_node -1;
                end
                final_index = strt + index -1;
                Ap(final_index,curr_node) = 1;
                %disp(Ap)
                if(rem(length(path),2)==0)
                    if path(j+2)==n
                        break
                    end
                    
                else
                    if next_node==n
                        break
                    end
                    
                end
            end
            flag = ~flag;
            
            
            
         end
    %end
    %fprintf("E is :%d\n",E)
    path_arr=reshape(Ap,[E,1]);
end