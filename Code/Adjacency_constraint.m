function opt = Adjacency_constraint(x_n_sq,B,C,h,src_id)

        A_mat = diag(x_n_sq);
        Z = B*A_mat*C;
        %display(Z);
        nn=[]; 
        for nd=1:h
            Za = Z^nd;
            nn=[nn Za(end,src_id)];
        end
        
        c = zeros(1,h);
        c(h)=1;
           
        if nn==c
            opt=1;
        else
            opt=0;
        end
        %end
end