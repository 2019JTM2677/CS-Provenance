function [index,g11,g12,g31,g32,reverse_group,edges] = max_index(A,Residue,g11,g12,g31,g32,reverse_group,edges,prev_index,n,alpha)
% max_index: Selects maximum index based on values on those indices
% The values are (projection values)*(path values)
% Input
%	A           Sensing matrix
%	Residue     Residue vector of length m, updated in each itertion
%   edges       Matrix with each row having [a b] for a->b or a->...->b
%   g11         Group-1: Matrix with each row containing indices for either
% incoming or outgoing edges from open ends of the path created by edges 
%   g12         Group-1: Array containing incomig edges to corresponding to
%   [a b] present in edges when b=n
%   g31         Group-3: Matrix with indices correspoding to edges to discard, if
%   present solution cannot be a path,i.e., a-> & ->b  if edges has [a b]
%   g32         Group-3: Only a->, as b=n
%   reverse_group Group-3 if a->b present, b->a discarded
%   prev_index  Index selected in previous iteration
%   alpha       Values assigned to Group-1=alpha and Group-2=1-alpha
% Group-2: indices unrelated to edges present
% Outputs
%	index        Index selected

    %proj_val = normalize(abs(A' * Residue),'range');
    proj_val = abs(A' * Residue);
    
%--------- Find node_a->node_b for index selected in previous iteration ---------------%
    if rem(prev_index,(n-1)) == 0    % node_b = dest
        node_a = prev_index/(n-1);
        node_b = n;
    else                             % node_b != dest
        node_b = ceil(prev_index/(n-1));
        r = rem(prev_index,(n-1));
        if r < node_b
            node_a = r;
        else
            node_a = r +1;
        end
    end
    
    % Find if the index selected is present in group 1 
    d1 = ismember(prev_index,g11);
    d2 = ismember(prev_index,g12);
    
    % Update the search groups only if d1 or d2 is true
    if d1
        %fprintf("G1\n");
        
        [row,col]=find(g11==prev_index) ;% But g11 has one less number of row than edges in total=>solved using 0000
        
        % Update edges and groups
        if length(row)==1
            if col<=n-2 % If last edge is [a_last b_last], then col<=n-2 is for ->a ie b_last=a, full edge= a_last->b
                new_a = node_a ;
                new_b = edges(row,2);
            else % for b-> then a_last = b, full edge= a->b_last, here b_last can be n
                new_a = edges(row,1);
                new_b = node_b; %edges(end,2);
            end
        elseif length(row)==2
            if (col(1)<=n-2 && col(2)>n-2) 
                new_a = edges(row(2),1);
                new_b = edges(row(1),2);
            elseif (col(2)<=n-2 && col(1)>n-2)
                new_a = edges(row(1),1);
                new_b = edges(row(2),2);
            end
        end
        edges(row,:)=[];
        try
            edges = cat(1,edges,[new_a new_b]);
        catch
           display(edges);
           g=[new_a new_b];
           display(g);
        end    
        %display(edges);
        g11(row,:)=[];
        %reverse_group(row,:) = [];
        %%g31(row,:)=[];
        
        if new_b == n
            
            % Final link new_a->new_b 
            g12 = cat(1,g12,Incoming_edges(new_a,n)); 
            %g12 = setdiff(g12,reverse_group);
            g12=g12 .* (~ismember(g12,reverse_group));
            
            
            % Zero for ->n link in g11 and g31 which contains ->!n
            g11 = cat(1,g11,zeros(1,(n-1)+(n-2)));
            reverse_group = cat(1,reverse_group,0);
            g31 = cat(1,g31,zeros(1,(n-1)+(n-2)));
            
            %--------- For group 3 : new_a-> & ->new_b, b->a --------------%
            search_spc=[];
            
            %corresponding to ->new_b
            search_spc = Incoming_edges(node_b,n);
            
            spc=[];
            
            %corresponding to new_a->
            spc = Outgoing_edges(node_a,n);
            search_spc = [search_spc spc];
            g32 = cat(1,g32,search_spc);
            % ---------------------------------------------------------------------%
        else
            %corresponding to ->new_a
            search_spc=[];
            search_spc = Incoming_edges(new_a,n);
            %corresponding to new_b->
            spc=[];
            spc = Outgoing_edges(new_b,n);
            search_spc = [search_spc spc];
            % find index for node_b->node_a  & new_b->new_a & remove from search space
            reverse_group = cat(1,reverse_group,reverse_edge(new_a,new_b,n));
            
            reverse_group = cat(1,reverse_group,reverse_edge(node_a,node_b,n));

            search_spc = search_spc .* (~ismember(search_spc,reverse_group));
            %display(search_spc)
            g11 = cat(1,g11,search_spc);

            %--------- For group 3 : a-> & ->b, b->a ---------------------%
            search_spc=[];
            %corresponding to ->new_b
            search_spc = Incoming_edges(node_b,n);
            %corresponding to new_a->
            spc=[];
            spc = Outgoing_edges(node_a,n);
            search_spc = [search_spc spc];
            g31 = cat(1,g31,search_spc);
            % ------------------------%
        
        end
        
    elseif d2  % this is related to a->n, so g contain ->a. So blast=a, hence final edge is alast->n
        
        %fprintf("G2\n");
        % [row,col]=find(g12==index);
        try
            [row,~] =find(g11==zeros(1,n-1+n-2)); % Find all zero row in g11 or just 0 in reverse group
        catch
            display(g11);
        end    
        new_a = node_a;%edges(end,1);
        new_b = n;
        %fprintf("Updates: %d %d\n",new_a,new_b);
        edges(row,:)=[];
        %%g31(row,:)=[];
        
        % Remove end row and all zero row from g11 and reverse
        g11(row,:)=[];
        reverse_group(reverse_group==0,:) = [];
        
        reverse_group = cat(1,reverse_group,reverse_edge(node_a,node_b,n));
        
        g12 = [];
        g32 = [];
        
        % Add zeros at last corresponding to new_a->n
        g11 = cat(1,g11,zeros(1,(n-1)+(n-2)));
        reverse_group = cat(1,reverse_group,0);
        g31 = cat(1,g31,zeros(1,(n-1)+(n-2)));
        
        edges = cat(1,edges,[new_a new_b]);
        
        search_spc = Incoming_edges(new_a,n);
        g12 = cat(1,g12,search_spc);
        
        %------------ For group 3 : a-> & ->b, b->a --------------------------%
        search_spc=[];

        %corresponding to ->new_b
        search_spc = Incoming_edges(new_b,n);

        %corresponding to new_a->
        spc=[];
        spc = Outgoing_edges(node_a,n);
        
        search_spc = [search_spc spc];
        g32 = cat(1,g32,search_spc);
        % ------------------------------------------------------------------%

        %}
    else
        %fprintf("Not match\n");
        new_a=node_a; new_b=node_b;
        edges = cat(1,edges,[new_a new_b]);
        if new_b == n
            g12 = cat(1,g12,Incoming_edges(new_a,n)); 
            g11 = cat(1,g11,zeros(1,(n-1)+(n-2)));
            reverse_group = cat(1,reverse_group,0);
            g31 = cat(1,g31,zeros(1,(n-1)+(n-2)));
            %------------ For group 3 : a-> & ->b, b->a --------------------------%
            search_spc=[];
            %corresponding to ->new_b
            search_spc = Incoming_edges(new_b,n);
            %corresponding to new_a->
            spc=[];
            spc = Outgoing_edges(new_a,n);
            
            search_spc = [search_spc spc];
            g32 = cat(1,g32,search_spc);
            % ---------------------------------------------------------------------%
        else
            %corresponding to ->new_a
            search_spc=[];
            search_spc = Incoming_edges(new_a,n);
            
            %corresponding to new_b->
            spc=[];
            spc = Outgoing_edges(new_b,n);
            search_spc = [search_spc spc];
            
            % Find index for node_b->node_a & remove from search space\
            
            ind = reverse_edge(new_a,new_b,n);

            search_spc(abs(search_spc)==ind)=0;%[];
            reverse_group = cat(1,reverse_group,ind);
            g11 = cat(1,g11,search_spc);

            %------------ For group 3 : a-> & ->b, b->a --------------------------%
            search_spc=[];
            %corresponding to ->new_b
            search_spc = Incoming_edges(node_b,n);%1+ (new_b - 1)*(n-1): new_b*(n-1) -1  ;
            %corresponding to new_a->
            spc=[];
            spc=Outgoing_edges(node_a,n);
            
            search_spc = [search_spc spc];
            g31 = cat(1,g31,search_spc);
            % ------------------------%
        
        end
    end
    
    %---------------- Create final array -------------------------------%
    % Assign values based on path likelihood
    path_val=[];
    for i=1:(n-1)^2
        if ismember(i,g31) || ismember(i,reverse_group)|| ismember(i,g32) % discarded edges
           path_val(i) = 0; 
        elseif ismember(i,g11) || ismember(i,g12)    % connected edges
           path_val(i) = alpha;
        else % disjoint edges
           path_val(i) = 1-alpha;
       end 
    end
    
    % Final values
    val = path_val' .* proj_val;
    
    [~,index] = sort(val,'descend');
    index = index(1);
    
   
end