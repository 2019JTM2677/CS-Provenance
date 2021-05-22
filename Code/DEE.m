% Double edge embedding
%function []= DEE(path,n,h)
% Functions used: node, path array de, double edge embed, OMP, OMP modified, 
%DE to SE, DE valid path, packet
    clear;
    n = 7; % no. of nodes
    hop_len = [1 2 3 4 5 6 7 8]; % Required hop length
    h = hop_len(5);
    if rem(h+1,2)==0    
        h_de = (h+1)/2 -1;  % Sparsity for double edge vector
    else
        h_de = floor((h+1)/2);
    end
    topology = 'complete';      % Or 'random'
    pkt_path = [2 3 4 1 5 n];
    no_of_pkts = 10^6;          % No. of iterations
    error_threshold = 200;      % Error Rate = error_threshold/no_of_pkts
    % Six values for column size 'm'
    start = 2*h_de;
    no_of_val = 8;
    inc = 5;  % Increment
    stop = min((no_of_val -1)*inc + start,((n-2)^2)*(n-1)); %m<DE
    m = [start:inc:stop];        % Column size values using AP
    N = min(min(m)/h_de,h_de);   % N for gOMP
    L = N+1;                     % List size L-OMP
    m=[45 50];
    
    % Run algorithm if value is set 1
    omp=0;
    gomp=0;
    cvx=0;
    Lomp=1;
    PLomp=0;
    Lgomp=1;
    % Save results in mat file if set to 1
    save_result=0;
    
    %% %--------------------- Form Topology -------------------------------%
    switch lower(topology)
    case 'complete'
        graph_adj_mat = ones(n) - diag(ones([1,n])); % adjacency matrix for complete graph
        connected_G = graph(graph_adj_mat~=0);       % form graph using adj mat
        path = pkt_path;                             % assign path
    case 'random'
        l = 100; b = 100; %length breadth
        points = [];
        for i=1:n
            points = [points; [randi(l,1,1) randi(b,1,1)]];
        end

        % Distance between coordinates: Adjacency Matrix
        Adj = squareform(pdist(points,'euclidean'));
        G = graph(Adj);

        [T,pred] = minspantree(G);       % Form MST using graph G
        maxweight = max(T.Edges.Weight); % Max edge weight of MST
        reachable_adj_mat = Adj; 
        
        % Form connected graph using max edge weight
        reachable_adj_mat(reachable_adj_mat <= maxweight+5 & reachable_adj_mat >0) = 1;
        reachable_adj_mat(reachable_adj_mat > maxweight+5) = 0;
        graph_adj_mat = reachable_adj_mat;
        %{ 
        % To fix a random topology
        graph_adj_mat =  [0     0     0     1     1     1     1;
         0     0     0     1     0     1     0;
         0     0     0     1     0     0     1;
         1     1     1     0     0     1     0;
         1     0     0     0     0     0     1;
         1     1     0     1     0     0     0;
         1     0     1     0     1     0     0];
        %}
        connected_G = graph(graph_adj_mat);
        
        % Select path of required hop length 'h'
        path=[];
        for i=1:n-1
            path = allpaths(connected_G,i,n,'MinPathLength',h,'MaxPathLength',h);
            if ~isempty(path)
                path = path{1};
                break
            end
        end
        
    end
    figure(1);
    plot(connected_G);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Simulation Parameters
    DE = ((n-2)^2)*(n-1); % No. of Double edges
    EE = n^2; 
    
    mu = 0; sigma = 1;  % gaussian parameters
    
    B=[];
    for ii=1:n
        count=0;
        for j=1:EE
            k=n*(ii-1)+1;
            count=k+n;
            if j>=k && j<count 
                B(ii,j)=1;
            else 
                B(ii,j)=0;
            end
        end
    end
    C=[];
    for ii=1:n
        C=[C;eye(n)];
    end

    temp_a=zeros(1,n^2);temp_b=zeros(1,n^2);
    temp_c=zeros(1,n^2);temp_d=zeros(1,n^2);
    ii=0;
    for j=1:n
        for i=1:n
            ii=ii+1;
            temp_a(ii)= n^2*(j-1)+ n*(j-1)+i;
            temp_b(ii)= n^2*j -i + 1;
            temp_c(ii)= n^2*(j-1)+ j + n*(i-1);
            temp_d(ii)= n^2*(j-1)+ (n+1)*(i-1)+1;

        end
    end
    temp_e = n^2*(n-1)+(1:n^2);
    A_1=zeros(n^2,n^3); %a->b->i
    A_2=zeros(n^2,n^3); %i->a->b

    for i=1:n^2
        j=(1:n)+n*(i-1);
        A_1(i,j)=1;
    end

    for i=1:n^2
        g=ceil(i/n);
        k= mod(i,n);
        if k==0
            k=n;
        end
        j=g+(0:n-1)*n+(n^2*(k-1));
        A_2(i,j)=1;
    end
            
    % Define n nodes
    nodes_de=[];  % All varE used for additional constraint
    for i= 1:n
        nodes_de = [nodes_de node(i,[],graph_adj_mat(i,:))];
    end

    % Destination
    dest = nodes_de(n);
    %fprintf("Destination node: %d",dest.Node_id);
    src = nodes_de(path(1));
    penultimate = nodes_de(path(h)); % Penultimate node
    
    fprintf('\nPath selected:');disp(path)
    
    % Path array similar to x (y=Ax) for verify
    Path_arr_de = path_array_de(path,DE,n);
    
    fprintf("Double edge embedding\n");
    
    error_rate=zeros(1,length(m));   %DE CVX
    error_rate_OMP=zeros(1,length(m));  % DE OMP
    error_rate_OMP_mod=zeros(1,length(m)); % DE modified OMP with brute force
    error_rate_OMP_mod_v1=zeros(1,length(m));  % Module 1 (missing double edge check)
    error_rate_OMP_mod_v2=zeros(1,length(m)); % Module 2   (Path constraints check)
    error_rate_OMP_topo = zeros(1,length(m)); % SE OMP with topology knowledge
    error_rate_gOMP = zeros(1,length(m));       % SE generalised OMP
    error_rate_gOMP_mod_v2 = zeros(1,length(m));
   
    for m_index =1: length(m)
        error_count=0;
        error_count_OMP=0;
        error_count_OMP_mod_v1=0;
        error_count_OMP_mod_v2=0;
        error_count_OMP_mod=0;
        error_count_OMP_topo =0;
        error_count_gOMP=0;
        error_count_gOMP_mod_v2=0;
        
        error_rate(m_index) = 0;
        error_rate_OMP(m_index) = 0;
        error_rate_OMP_mod_v1(m_index) = 0;
        error_rate_OMP_mod_v2(m_index) = 0;
        error_rate_OMP_mod(m_index) = 0;
        error_rate_OMP_topo(m_index)=0;
        error_rate_gOMP(m_index)=0;
        error_rate_gOMP_mod_v2(m_index)=0;
        
        pkt_count = 0;

        %sigma=1/m(m_index);
        for pkt_i=1:no_of_pkts
            
            pkt =  Packet;
            pkt_count = pkt_count +1;
            Ar = normrnd(mu,sigma,[m(m_index),DE]);
            
            % Assigning columns to edges
            Edge_id_de = reshape(Ar,[m(m_index),(n-2)^2,n-1]);

            for ni= 1:n
                if ni==n
                    dbl_edge = Ar;
                else
                    dbl_edge = Edge_id_de(:,:,ni);
                end
                nodes_de(ni).Edge_id = dbl_edge;
            end

            % Embed double edge ids
            pkt = double_edge_embed(m(m_index),path,pkt,nodes_de);
            y = pkt.provenance;
            
%-------------------- Recovery using OMP double edge ---------------------%
            if omp
                x_OMP = OMP(h_de,y,Ar);
                x_OMP(abs(x_OMP)<=0.001)=0;
                x_OMP(abs(x_OMP)>0.001)=1;

                if ~isequal(x_OMP, Path_arr_de)
                   error_count_OMP = error_count_OMP +1; % Increment count
                end
            end
%-------------------------------------------------------------------------%            
            
%---------------------- Recovery using gOMP ------------------------------%
            if gomp
                %N = 4;%floor(min(h_de,m(m_index)/h_de));
                x_gOMP = gOMP(h_de,y,Ar,N);
                x_gOMP(abs(x_gOMP)<=0.001)=0;
                x_gOMP(abs(x_gOMP)>0.001)=1;
                if ~isequal(x_gOMP, Path_arr_de)
                   error_count_gOMP = error_count_gOMP +1; % Increment count
                end
            end
%-------------------------------------------------------------------------%            
            

%---------------------- Recovery using modified OMP module 2 -------------%
            %{
            % Verifying using adajacency matrix
            B=[];
            for ii=1:n
                count=0;
                for j=1:EE
                    k=n*(ii-1)+1;
                    count=k+n;
                    if j>=k && j<count 
                        B(ii,j)=1;
                    else 
                        B(ii,j)=0;
                    end
                end
            end
            C=[];
            for ii=1:n
                C=[C;eye(n)];
            end

            temp_a=zeros(1,n^2);temp_b=zeros(1,n^2);
            temp_c=zeros(1,n^2);temp_d=zeros(1,n^2);
            ii=0;
            for j=1:n
                for i=1:n
                    ii=ii+1;
                    temp_a(ii)= n^2*(j-1)+ n*(j-1)+i;
                    temp_b(ii)= n^2*j -i + 1;
                    temp_c(ii)= n^2*(j-1)+ j + n*(i-1);
                    temp_d(ii)= n^2*(j-1)+ (n+1)*(i-1)+1;

                end
            end
            temp_e = n^2*(n-1)+(1:n^2);
            A_1=zeros(n^2,n^3); %a->b->i
            A_2=zeros(n^2,n^3); %i->a->b

            for i=1:n^2
                j=(1:n)+n*(i-1);
                A_1(i,j)=1;
            end

            for i=1:n^2
                g=ceil(i/n);
                k= mod(i,n);
                if k==0
                    k=n;
                end
                j=g+(0:n-1)*n+(n^2*(k-1));
                A_2(i,j)=1;
            end
            %}
            if Lomp
               
                OMP_ts_module2 = tic;
                x_OMP_mod_v2 = OMP_modified(h_de,y,Ar,L);
                x_OMP_mod_v2(abs(x_OMP_mod_v2)<=0.001)=0;
                x_OMP_mod_v2(abs(x_OMP_mod_v2)>0.001)=1;

                size_OMP_v2 = size(x_OMP_mod_v2,2);

                jj=[];
                for j=1:size(x_OMP_mod_v2,2) 
                    x_n_sq = DE_to_SE(x_OMP_mod_v2(:,j),n,h,penultimate.Node_id,temp_a,temp_b,temp_c,temp_d,temp_e,A_1,A_2);
                    opt = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
                    if opt==0
                            jj=[jj j];
                    end
                end
                x_OMP_mod_v2(:,jj)=[];
                OMP_te_module2(pkt_i) = toc(OMP_ts_module2);

                %Arrange in ascending residual norm
                residual_norm=zeros(1,size(x_OMP_mod_v2,2));
                for j=1:size(x_OMP_mod_v2,2)
                    final_residue = y - Ar*x_OMP_mod_v2(:,j);
                    residual_norm(j) = norm(final_residue);
                end
                [~,order]= sort(residual_norm);
                x_OMP_mod_v2 = x_OMP_mod_v2(:,order);

                if isempty(x_OMP_mod_v2)
                    rec_x_OMP_mod_v2 = zeros(DE,1);
                else
                    rec_x_OMP_mod_v2 = x_OMP_mod_v2(:,1);
                end

                if ~isequal(rec_x_OMP_mod_v2, Path_arr_de)
                   error_count_OMP_mod_v2 = error_count_OMP_mod_v2 +1; % Increment count
                   flag_v2=0;
                   if ismember(Path_arr_de',x_OMP_mod_v2','rows')
                        fprintf("RESIDUE BUG\n") 
                   end
                else
                    flag_v2=1;
                end
            end
%--------------- verifying using brute force ------------------------------%
            %{
            ts_brute = tic;
            ji=[];
            for j=1:size(x_OMP_mod,2)
                opt = DE_valid_path(x_OMP_mod(:,j),n,h,src.Node_id,dest.Node_id,penultimate.Node_id);
                if opt==0
                    ji=[ji j];
                end
            end
            x_OMP_mod(:,ji)=[];
            te_brute(pkt_i) = toc(ts_brute) ;

            if isempty(x_OMP_mod)
                rec_x_OMP_valid = zeros(DE,1);%orig_x_OMP_mod(:,1);
            else
                rec_x_OMP_valid = x_OMP_mod(:,1);
            end

            if ~isequal(rec_x_OMP_valid, Path_arr_de)
               error_count_OMP_mod = error_count_OMP_mod +1; % Increment count
            end
            
            if ~isequal(rec_x_OMP_adj,Path_arr_de) && isequal(rec_x_OMP_valid,Path_arr_de)
                fprintf("bug");
                break
            end
            %}
%-------------------------------------------------------------------------%

            %
%----------------------- modified OMP module 1 ---------------------------%
            if PLomp
                
                OMP_ts_module1 = tic;
                x_OMP_mod_v1 = PL_OMP_DE(n,h_de,h,y,Ar,L,B,C,A_1,A_2,temp_a,...
                temp_b,temp_c,temp_d,temp_e,src.Node_id,dest.Node_id,penultimate.Node_id);    
                OMP_te_module1(pkt_i) = toc(OMP_ts_module1);

                size_OMP_v1 = size(x_OMP_mod_v1,2);
                
                if isempty(x_OMP_mod_v1)
                    rec_x_OMP_mod_v1 = zeros(DE,1);
                else
                    rec_x_OMP_mod_v1 = x_OMP_mod_v1(:,1);
                end

                if ~isequal(rec_x_OMP_mod_v1, Path_arr_de)
                    error_count_OMP_mod_v1 = error_count_OMP_mod_v1 +1; % Increment count
                end
            end
%-------------------------------------------------------------------------%
            %}
%-------------------- Recovery using L-gOMP -----------------------%
            if Lgomp
                %N = 4;%floor(min(h_de,m(m_index)/h_de));
                x_gOMP_mod_v2 = gOMP_list_module2(h_de,y,Ar,N);

                x_gOMP_mod_v2(abs(x_gOMP_mod_v2)<=0.001)=0;
                x_gOMP_mod_v2(abs(x_gOMP_mod_v2)>0.001)=1;
                orig_x_gOMP_mod=x_gOMP_mod_v2;
                size_gOMP_v2 = size(x_gOMP_mod_v2,2);

                jj=[];              % Invalid path indices from set of solutions
                for j=1:size(x_gOMP_mod_v2,2)
                        x_n_sq = DE_to_SE(x_gOMP_mod_v2(:,j),n,h,penultimate.Node_id,temp_a,temp_b,temp_c,temp_d,temp_e,A_1,A_2);
                        opt = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
                        if opt==0
                            jj=[jj j];
                        end
                end
                x_gOMP_mod_v2(:,jj)=[];

                % Arrange in ascending residual norm
                residual_norm=zeros(1,size(x_gOMP_mod_v2,2));
                for j=1:size(x_gOMP_mod_v2,2)
                    final_residue = y - Ar*x_gOMP_mod_v2(:,j);
                    residual_norm(j) = norm(final_residue);
                end
                [~,order]= sort(residual_norm);
                x_gOMP_mod_v2 = x_gOMP_mod_v2(:,order);

                if isempty(x_gOMP_mod_v2)
                    rec_x_gOMP_mod_v2 = zeros(DE,1);%orig_x_gOMP_mod(:,1);
                else
                    rec_x_gOMP_mod_v2 = x_gOMP_mod_v2(:,1);
                end

                if ~isequal(rec_x_gOMP_mod_v2, Path_arr_de)
                   error_count_gOMP_mod_v2 = error_count_gOMP_mod_v2 +1; % Increment count
                   if ismember(Path_arr_de',orig_x_gOMP_mod','rows')
                        fprintf("RESIDUE BUG\n") 
                   end
                end
            end
%-------------------------------------------------------------------------%            
            %}
            
            %{
%------------- Recovery using OMP when topology is known -----------------%
            % Selecting columns of Ar using topology adjacency matrix
            temp_k=[];
            kk=1;
            for i=1:(n-1)
                for j=1:(n-1)
                    if i==j
                    elseif graph_adj_mat(i,j)==1
                        for b=1:n
                            if (b==i) || (b==j)
                                %fprintf("Not embedding:%d %d %d %d",i,j,b,graph_adj_mat(b,i));
                            elseif graph_adj_mat(b,i)==1
                                temp_k=[temp_k kk];
                                kk=kk+1;
                                %fprintf("Embed: %d %d %d %d",i,j,b,graph_adj_mat(b,i));
                            else
                                kk=kk+1;
                            end
                        end
                    else 
                        for l=1:n-2
                            kk=kk+1;
                        end
                    end
                end
            end
            A_topo = Ar(:,temp_k);
            x_topo = OMP_modified(h_de,y,A_topo,L);
            x_topo(abs(x_topo)<=0.001)=0;
            x_topo(abs(x_topo)>0.001)=1;
            
            x_OMP_top=zeros(DE,size(x_topo,2));
            for uu=1:size(x_topo,2)
                x_OMP_top(temp_k,uu) = x_topo(:,uu);
            end
            
            orig_x_OMP_top = x_OMP_top;
            ji=[];
            for j=1:size(x_OMP_top,2) 
                x_n_sq = DE_to_SE(x_OMP_top(:,j),n,h,penultimate.Node_id,temp_a,temp_b,temp_c,temp_d,temp_e,A_1,A_2);
                opt1 = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
                
                A_mat = diag(x_n_sq);
                Z = B*A_mat*C;
                nn=zeros(1,h); 
                for hi=1:h
                    Za = Z^hi;
                    nn(hi)=Za(end,src.Node_id);
                end
                c = zeros(1,h);
                c(h)=1;
                if nn==c 
                   %fprintf("Valid h len path");
                   if opt1~=1
                       fprintf("BUG\n");
                       return
                   end
                else
                    %fprintf("Invalid");
                    ji=[ji j];
                end
            end
            x_OMP_top(:,ji)=[];

            if isempty(x_OMP_top)
                rec_x_OMP_top = zeros(DE,1);%orig_x_OMP_top(:,1);
            else
                rec_x_OMP_top = x_OMP_top(:,1);
            end
            
            if ~isequal(rec_x_OMP_top, Path_arr_de)
               error_count_OMP_topo = error_count_OMP_topo +1; % Increment count
            end
%-------------------------------------------------------------------------%
            %}
%------------------- Recovery using CVX double edge ----------------------%
            
            %x_de = cvx_solver(DE,b_de,Ar_de);
            if cvx
                cvx_ts_module = tic;
                cvx_begin quiet
                cvx_precision default
                variable x_de(DE) %binary
                minimize( norm(x_de, 1 ) )
                subject to
                    Ar*x_de == y;
                cvx_end

                x_de(abs(x_de)<=0.001)=0;
                x_de(abs(x_de)>0.001)=1;
                cvx_te_module(pkt_i) = toc(cvx_ts_module) ;

                if ~isequal(x_de,Path_arr_de)
                    error_count = error_count +1; %Increment if path recovered is different from path travelled
                end
            end    
          
%-------------------------------------------------------------------------%
            % HALT CONDITION
            if error_count_gOMP_mod_v2 == error_threshold
                break
            end

        end
        
        % Time
        %{
        OMP_t_module1(m_index) = mean(OMP_te_module1);
        OMP_t_module2(m_index) = mean(OMP_te_module2);
        cvx_t_module(m_index) = mean(cvx_te_module);
        fprintf("Time complexity:") ;
        display(OMP_t_module1);
        display(OMP_t_module2);
        %}
        
        % List size
        %fprintf("size list OMP module 1:");display(size(x_OMP_mod_v1));
        %fprintf("\nsize list L-OMP:%d\n",size_OMP_v2);
        %fprintf("size list L-gOMP:%d\n",size_gOMP_v2);       
        
        error_rate(m_index) = error_count/pkt_count;
        error_rate_OMP(m_index) = error_count_OMP/pkt_count;
        error_rate_OMP_mod_v1(m_index) = error_count_OMP_mod_v1/pkt_count;
        error_rate_OMP_mod_v2(m_index) =  error_count_OMP_mod_v2/pkt_count;
        %error_rate_OMP_mod(m_index) =  error_count_OMP_mod/pkt_count;
        %error_rate_OMP_topo(m_index) = error_count_OMP_topo/pkt_count;
        error_rate_gOMP(m_index) =  error_count_gOMP/pkt_count;
        error_rate_gOMP_mod_v2(m_index) =  error_count_gOMP_mod_v2/pkt_count;
        
        %fprintf("Error Rate:%f for column size %d\n",error_rate(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size OMP %d\n",error_rate_OMP(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size PL-OMP %d\n",error_rate_OMP_mod_v1(m_index),m(m_index));
        fprintf("Error Rate:%f for column size L-OMP %d\n",error_rate_OMP_mod_v2(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size mod OMP %d\n",error_rate_OMP_mod(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size OMP topo %d\n",error_rate_OMP_topo(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size gOMP  %d\n",error_rate_gOMP(m_index),m(m_index));
        fprintf("Error Rate:%f for column size L-gOMP %d\n\n",error_rate_gOMP_mod_v2(m_index),m(m_index));
        
    end
    % end

    %% PLOTS
    %
    figure(2)
    semilogy(m,error_rate, 'yo-', 'LineWidth', 1.5);
    hold on
    semilogy(m,error_rate_OMP, 'b+-', 'LineWidth', 1.5);
    semilogy(m,error_rate_OMP_mod_v1, 'ms-', 'LineWidth', 1.5);
    semilogy(m,error_rate_OMP_mod_v2, 'gs-', 'LineWidth', 1.5);
    %semilogy(m,error_rate_OMP_topo, 'r+-', 'LineWidth', 1.5);
    semilogy(m,error_rate_gOMP, 'm+-', 'LineWidth', 1.5);
    semilogy(m,error_rate_gOMP_mod_v2, 'ks-', 'LineWidth', 1.5);
    axis([0 120 0 1]);
    grid on
    legend('CVX','OMP','L-OMP','gOMP','L-gOMP');
    title('Error rate vs column size DE');
    xlabel('Column size');
    ylabel('Error rate');
    hold off
%end
%}
infos ={'Results for SEE:(I)CVX (II)1. OMP 2. L-OMP (III)1. gOMP 2. L-gOMP';
   'n= no. of nodes';'h = path hop length';'h_de = DE sparsity';'DE= no. of possible double edges';
   'm = column size of Ar for OMP';'graph_adj_mat = topology adjacency matrix'};
if save_result
   save DEE_n9_h8.mat infos n DE h h_de graph_adj_mat path m N error_rate error_rate_OMP error_rate_OMP_mod_v2 error_rate_gOMP error_rate_gOMP_mod_v2
end