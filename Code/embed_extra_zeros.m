function [x_n_sq]= embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest)

    u=1;v=1;
    x_i2_dest = x_nminus1_sq(temp_dest);
    x_nminus1_sq(temp_dest)=[];
    x_n_sq = zeros(n^2,1);
    for k=1:n^2
        if (ismember(k,temp_e)) || (ismember(k,temp_f))
             x_n_sq(k)=0;    
        elseif (ismember(k,temp_g))
             x_n_sq(k)= x_i2_dest(u);
             u=u+1;
        else
             %fprintf("k=%d",k);
             %fprintf(" v=%d\n",v);
             x_n_sq(k)= x_nminus1_sq(v);
             %v = mod(v+1,length(x_nminus1_sq)+2);
             v=v+1;
        end
        %fprintf("k=%d",k);
    end

end