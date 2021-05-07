% Function to identify recovered solution is valid path,opt=1 for valid
function opt = valid_path(x,n,h,src_id)
       
       index = find(x==1);     % non zero indices
       if size(index,1)~=h
           opt=0;
       else    
           dest_ind=zeros(1,n-1);
           for i=1:(n-1)
               dest_ind(i)= i*(n-1);  % indices for a->b, b=dest
           end

           opt=1;
           nonzero_dest = intersect(index,dest_ind);
           c=[];
           if length(nonzero_dest)==1     % exactly 1 non zero index for i->dest
               for i = 1: length(index)
                    if mod(index(i),n-1)==0
                        a = index(i)/(n-1) ;  %a->b
                        b = n;
                        c=[c;[a,b]];
                    else
                        b1 = ceil(index(i)/(n-1));
                        a1 = mod(index(i),n-1);
                        if a1 >= b1
                            a1=a1+1;
                        end
                        c=[c;[a1,b1]];

                    end
                end
                bd = c(:,2); 
                ad = c(:,1);
                cd = c;
                ii=1;
                path=[];
                if length(bd) ~= length(unique(bd)) || length(ad) ~= length(unique(ad))
                    opt=0;

                else
                    while ~isempty(cd)
                        if bd(ii)==b
                            if size(cd,1)==1
                                if ad(ii)==src_id  %src.Node_id=1
                                    opt=1;
                                    path=[src_id b path];
                                else
                                    opt=0;
                                end
                            else    
                                opt=1;
                                path=[b path];
                            end
                            %fprintf("found %d\n",ii)
                            b = ad(ii);
                            cd(ii,:)=[];
                            %fprintf("size:%d\n",size(cd,1));
                            bd = cd(:,2) ;
                            ad = cd(:,1);
                            ii=1;

                        else
                            opt = 0;
                            %fprintf("Not found %d\n",ii)
                            ii=ii+1;

                            if ii==size(cd,1)+1
                                %fprintf("size %d\n",size(cd,1));
                                break
                            end
                        end
                    end
                end
                %display(path);
            else % no index corresponding to ->dest =1
               opt = 0;
           end
       end
end