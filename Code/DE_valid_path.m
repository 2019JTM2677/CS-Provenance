% Function to identify recovered solution is valid path,opt=1 for valid
function opt = DE_valid_path(x,n,h,src_id, dest_id, penult_id)
       
       index = find(x==1);     % non zero indices
       path= [];%zeros(1,h+1);
       
       if size(index,1)~=floor(h/2)
           opt=0;
          
       else    
           dest_ind=zeros(1,(n-2)*(n-1));
           for i=1:(n-1)*(n-2)
               dest_ind(i)=i*(n-2);  % indices for a->b, b=dest
           end
           
           opt=1;
           nonzero_dest = intersect(index,dest_ind);
           d=[];
           if (length(nonzero_dest)==1  && mod(h,2)==0) || mod(h,2)==1   % exactly 1 non zero index for i->dest
               %Finding a,b,c from non zero index position
               for i = 1: length(index)
                    b = ceil(index(i)/(n-2)^2); %a->b->c
                    
                    x = ceil((index(i)-(b-1)*(n-2)^2)/(n-2));
                    if x >= b
                        a=x+1;
                    else
                        a=x;
                    end
                    
                    y = index(i) - (b-1)*(n-2)^2 -(n-2)*(x-1);
                    %display(y<a);
                    %display(y<b);
                    if y<a && y<b
                        c=y; 
                    elseif y+2 >a && y+2 >b
                        c=y+2;
                    else
                        c=y+1;
                    end
                    d=[d;[a,b,c]];

                end
                ad = d(:,1); 
                bd = d(:,2);
                cd = d(:,3);
                dd = d;
                %display(dd);
                ii=1;
                if mod(h,2)==0
                    temp_c = dest_id;
                else
                    temp_c = dd(end,3);
                end
                path=[];
                %Repeated node condition
                rep1 = ~isempty(intersect(ad,bd));
                rep2 = ~isempty(intersect(cd,bd));
                rep3 =  length(bd) ~= length(unique(bd)) || length(ad) ~= length(unique(ad)) || length(cd) ~= length(unique(cd)) ;
                if rep1 || rep2 || rep3 
                    opt=0;

                else
                    while ~isempty(dd)
                        if cd(ii)== temp_c 
                            if size(dd,1)==1
                                if ad(ii)==src_id  %src.Node_id=1
                                    opt=1;
                                    path=[src_id bd(ii) temp_c path];
                                else
                                    opt=0;
                                    path=[ad(ii) bd(ii) temp_c path];
                                end
                            else    
                                opt=1;
                                path=[bd(ii) temp_c path];
                            end
                            %fprintf("found %d\n",ii)
                            temp_c = ad(ii);
                            dd(ii,:)=[];
                            
                            cd = dd(:,3);
                            bd = dd(:,2) ;
                            ad = dd(:,1);
                            ii=1;

                        else
                            opt = 0;
                            %fprintf("Not found %d\n",ii)
                            ii=ii+1;

                            if ii==size(dd,1)+1
                                %fprintf("size %d\n",size(cd,1));
                                break
                            end
                        end
                    end
                end
                %display(path);
                % Odd hop length,last edge
                if mod(h,2)==1 
                    if path(end)== penult_id 
                        path=[path dest_id];
                    else
                        opt=0;
                    end
                end
                %{
                if isempty(path)
                    
                elseif path(end)~= dest_id %|| size(path,2)==h
                    path=[path dest_id];
                end
                %}
            else % no index corresponding to ->dest =1
               opt = 0;
           end
       end
       if isempty(path) || size(path,2)==1
           opt=0;
           
       elseif size(path,2)==h+1
           path_matrix = all_possible_paths(n,h,src_id,dest_id);
           %display(path);
           %display(size(path));
           %display(size(path_matrix));
           if ismember(path,path_matrix,'rows')
               opt=1;
               %fprintf("yes")
           else
               opt=0;
               %fprintf("no")
           end
       end 
end