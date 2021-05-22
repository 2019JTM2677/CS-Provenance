% Single edge embedding on a network of n nodes and h hops
%function []= SEE(n,m,h,no_of_pkts, graph)

% Functions used: node, path array, edge embed, OMP, OMP modified, valid path, packet
    %clc;
    clear;
    
    % Input Values
    n = 7;                      % No. of nodes
    hop_len = [1 2 3 4 5 6 7 8]; 
    h = hop_len(6);             % Required path length
    topology = 'complete';      % Or 'random'
    pkt_path = [1 3 2 4 5 6 n];
    no_of_pkts = 10^3;          % No. of iterations
    error_threshold = 100;      % Error Rate = error_threshold/no_of_pkts
    % Six values for column size 'm' starting at 2h
    start = 2*h;
    no_of_val = 6;
    inc = 4;  % Increment
    stop = min((no_of_val -1)*inc + start,(n-1)^2);%m<E
    m = [start:inc:stop];       % Column size values using AP
    N = min(min(m)/h,h) ;       % N for gOMP 
    L = N+1;                    % List size L-OMP
    
    % Run algorithm if value is set 1
    omp=1;
    gomp=1;
    stomp=0;
    cosamp=0;
    cvx=1;
    Lomp=1;
    PLomp=1;
    PLomp2=1;
    Lgomp=1;
    PLgomp=1;
    % Save results in mat file if set to 1
    save_result=0;
    
    %% %--------------------- Form Topology ------------------------------%
    
    switch lower(topology)
    case 'complete'
        graph_adj_mat = ones(n) - diag(ones([1,n])); % adjacency matrix for complete graph
        connected_G = graph(graph_adj_mat~=0);       % form graph using adj mat
        path = pkt_path;                             % assign path
    case 'random'
        % Node distribution
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
    %% %------------------- Simulation Parameters-------------------------%
    
    E = (n-1)^2; % Total no of directed edges
    EE = n^2;    % Modified no of directed edges for additional constraint
    
    mu=0;sigma=1;          % gaussian parameters
    
    % verifying using adjacency, convert (n-1)^2 -> n^2
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

    % Extra indices in n^2 vector
    temp_e=zeros(1,n);
    temp_f=zeros(1,n-1);
    temp_g=zeros(1,n-1);

    for i=1:n
        temp_e(i)=(i-1)*n + i;
    end

    for i=1:n-1
        temp_f(i)=i*n;
        temp_g(i)= n*(n-1)+i;
    end

    temp_dest=zeros(1,n-1); % indices mapped to i->dest
    for i=1:(n-1)
        temp_dest(i)=i*(n-1);
    end
    
    %{
    temp_m = 2*h; 
    m=[]; % form array for various values of no of rows of matrix Ar(m,n) 
    for ii=1:5
        m =[m temp_m+(ii-2)*temp_m/2];
    end
    %}
    
    % Define n nodes 
    nodes_1 =[]; 
    for i= 1:n
        nodes_1 =[nodes_1 node(i,[],graph_adj_mat(i,:))];
    end 
    
    % Destination
    dest = nodes_1(n);
    %fprintf("Destination node: %d",dest.Node_id);
    src = nodes_1(path(1));
    
    fprintf('\nPath selected:');disp(path) 
    
    % Path vector similar to x (y=Ax) to verify
    Path_arr = path_array(path,E,n); 
    
    fprintf("Single edge embedding\n");
    
    error_rate = zeros(1,length(m));            % CVX
    error_rate_OMP = zeros(1,length(m));        % OMP
    error_rate_OMP_mod_v1=zeros(1,length(m));   % PL-OMP one missing link(in-built path constraints)
    error_rate_OMP_mod_v2 = zeros(1,length(m)); % L-OMP(path constraints after solution)
    error_rate_OMP_mod = zeros(1,length(m));    % L-OMP  with brute force
    error_rate_OMP_topo = zeros(1,length(m));   % L-OMP with topology knowledge
    error_rate_gOMP = zeros(1,length(m));       % gOMP
    error_rate_gOMP_mod_v1=zeros(1,length(m));  % PL-gOMP 
    error_rate_gOMP_mod_v2 = zeros(1,length(m));% L-gOMP 
    error_rate_cosamp=zeros(1,length(m));       % CoSaMP
    error_rate_stOMP=zeros(1,length(m));        % StOMP
    list_wo_path_rate=zeros(1,length(m));       % OMP solution list have no valid paths
    error_rate_PLOMP2=zeros(1,length(m));       % PL-OMP with two missing link
    
    for m_index =1: length(m)
        error_count = 0;    
        error_count_OMP = 0;  
        error_count_OMP_mod_v1 = 0;
        error_count_PLOMP2=0;
        error_count_OMP_mod_v2 = 0;
        error_count_OMP_mod = 0; 
        error_count_OMP_topo = 0;
        error_count_gOMP = 0;
        error_count_gOMP_mod_v1 = 0;
        error_count_gOMP_mod_v2 = 0;
        error_count_cosamp = 0;
        error_count_stOMP = 0;
        list_wo_path_count = 0;
        
        error_rate(m_index) = 0;
        error_rate_OMP(m_index) = 0;
        error_rate_OMP_mod_v1(m_index) = 0;
        error_rate_OMP_mod_v2(m_index) = 0;
        error_rate_OMP_mod(m_index) = 0;
        error_rate_OMP_topo(m_index)=0;
        error_rate_gOMP(m_index)=0;
        error_rate_gOMP_mod_v1(m_index)=0;
        error_rate_gOMP_mod_v2(m_index)=0;
        error_rate_cosamp(m_index)=0;
        error_rate_stOMP(m_index)=0;
        error_rate_PLOMP2(m_index)=0;
        list_wo_path_rate(m_index)=0;
        
        pkt_count = 0;

        %sigma=1/m(m_index);  % gOMP paper theorem
        for pkt_i=1:no_of_pkts

            pkt =  Packet;
            pkt_count = pkt_count +1;
            Ar = normrnd(mu,sigma,[m(m_index),E]); 
            

            % Assigning columns to edges
            Edge_id = reshape(Ar,[m(m_index),E/(n-1),n-1]);

            % Distributing edge ids to nodes
            for ni= 1:n
                if ni==n
                    edge=Ar;    
                else
                    edge=Edge_id(:,:,ni);
                end
                nodes_1(ni).Edge_id = edge;
            end

            % Embed edge ids for path in packet
            pkt = edge_embed(m(m_index),path,pkt,nodes_1); 
            y = pkt.provenance;
            %fprintf("Final provenance\n");disp(y)

            
%-------------------- Recovery using OMP ---------------------------------%
            if omp
                x_OMP = OMP(h,y,Ar);
                x_OMP(abs(x_OMP)<=0.001)=0;
                x_OMP(abs(x_OMP)>0.001)=1;
              
                if ~isequal(x_OMP, Path_arr)
                   error_count_OMP = error_count_OMP +1; % Increment count
                end
            end
%-------------------------------------------------------------------------%
            
%-------------------- Recovery using gOMP --------------------------------%
            if gomp
                %N = 2;%floor(min(h,m(m_index)/h));
                x_gOMP = gOMP(h,y,Ar,N);
                x_gOMP(abs(x_gOMP)<=0.001)=0;
                x_gOMP(abs(x_gOMP)>0.001)=1;

                if ~isequal(x_gOMP, Path_arr)
                   error_count_gOMP = error_count_gOMP +1; % Increment count
                end
            end
%-------------------------------------------------------------------------%    
            
            %
%-------------------- Recovery using StOMP -------------------------------%
            if stomp
                %{
                delta = m(m_index)/E;
                rho = h/m(m_index);
                S = 10;
                alpha_0 = delta*(1-rho)/S;
                %q = min((m(m_index)-h)/h,0.5);
                [x_stOMP, iters] = SolveStOMP(Ar, y, E, 'FAR', alpha_0, S, 1);
                %}
                x_stOMP = StOMP(h,y,Ar);
                x_stOMP(abs(x_stOMP)<=0.001)=0;
                x_stOMP(abs(x_stOMP)>0.001)=1;
                if ~isequal(x_stOMP, Path_arr)
                   error_count_stOMP = error_count_stOMP +1; % Increment count
                else
                %     fprintf("bug")
                end
            end
%-------------------------------------------------------------------------%

%-------------------- Recovery using CoSaMP-------------------------------%
            if cosamp
                x_cosamp = CoSaMP(h,y,Ar);
                x_cosamp(abs(x_cosamp)<=0.001)=0;
                x_cosamp(abs(x_cosamp)>0.001)=1;
                if ~isequal(x_cosamp, Path_arr)
                   error_count_cosamp = error_count_cosamp +1; % Increment count
                end
            end
%-------------------------------------------------------------------------%            
            %}
            
            %}
%------------------- Recovery using OMP module 2-------------------------%
            %N=2;%floor(min(h,m(m_index)/h));
            %L=N+1;
            
            if Lomp
                OMP_ts_module2 = tic;       % Time elapse check
                [x_OMP_mod_v2,size_v2] = module_2(n,h,y,Ar,L,B,C,temp_e,temp_f,temp_g,temp_dest,src.Node_id);
                OMP_te_module2(pkt_i) = toc(OMP_ts_module2) ;

                if isempty(x_OMP_mod_v2)
                    rec_x_OMP_mod_v2 = zeros(E,1);
                    %list_wo_path_count = list_wo_path_count +1;
                else
                    rec_x_OMP_mod_v2  = x_OMP_mod_v2(:,1);
                end

                if ~isequal(rec_x_OMP_mod_v2, Path_arr)
                    %flag_v2=1;
                    error_count_OMP_mod_v2 = error_count_OMP_mod_v2 +1; % Increment count
                    if ismember(Path_arr',x_OMP_mod_v2','rows')
                        fprintf("RESIDUE BUG OMP\n") 
                   end
                else
                    %flag_v2=0;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            ts_module3 = tic;
            x_OMP_mod = OMP_modified2(h,y,Ar,L);
            x_OMP_mod(abs(x_OMP_mod)<=0.001)=0;
            x_OMP_mod(abs(x_OMP_mod)>0.001)=1;
            x_OMP_mod_adj = x_OMP_mod;
            %orig_x_OMP_mod = x_OMP_mod;
            
            jj=[];              % Invalid path indices from set of solutions
            for j=1:size(x_OMP_mod_adj,2)
                x_nminus1_sq = x_OMP_mod_adj(:,j);
                x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
                
                opt = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
                if opt==1
                   %fprintf("Valid h len path");
                else
                    %fprintf("Invalid");
                    jj=[jj j];
                end
                
            end
            x_OMP_mod_adj(:,jj)=[];
            te_module3(pkt_i) = toc(ts_module3) ;
            
            if isempty(x_OMP_mod_adj)
                rec_x_OMP_adj = zeros(E,1);%orig_x_OMP_mod(:,1);
                %list_wo_path_count = list_wo_path_count +1;
            else
                rec_x_OMP_adj = x_OMP_mod_adj(:,1);
            end
            %{
            % Bug check
            if ~isequal(x_OMP_mod_adj,x_OMP_mod_v2)
                fprintf("BUG")
                return
            end
            %}
            
            %{
            % verifying using brute force
            ji=[];
            %ts_brute = tic;
            for j=1:size(x_OMP_mod,2)
                opt = valid_path(x_OMP_mod(:,j),n,h,src.Node_id);
                if opt==1
                else
                    ji=[ji j];
                end
            end
            x_OMP_mod(:,ji)=[];
            
            % Bug check
            if ~isequal(x_OMP_mod,x_OMP_mod_v2)
                fprintf("BUG\n")
                return
            end
            
            %{
            te_brute(pkt_i) = toc(ts_brute) ;
            
            % Bug check
            if x_OMP_mod_adj ~= x_OMP_mod
                display(x_OMP_mod_adj);
                display(x_OMP_mod)
            end
            
            %}
            if isempty(x_OMP_mod)
                rec_x_OMP_valid = zeros(E,1);%orig_x_OMP_mod(:,1);
            else
                rec_x_OMP_valid = x_OMP_mod(:,1);
            end 
          %} 
            
            if ~isequal(rec_x_OMP_adj , Path_arr)
               error_count_OMP_mod = error_count_OMP_mod +1; % Increment count
            end
            %}
%-------------------------------------------------------------------------%    

%------------ Recovery using modified OMP module 1 (PL-OMP one missing link)----------%
            if PLomp
                OMP_ts_module1 = tic;
                x_OMP_mod_v1 = module_1(n,h,y,Ar,L,B,C,temp_e,temp_f,temp_g,temp_dest,src.Node_id,dest.Node_id,'SE');
                OMP_te_module1(pkt_i) = toc(OMP_ts_module1) ;

                if isempty(x_OMP_mod_v1)
                    rec_x_OMP_mod_v1 = zeros(E,1);
                else
                    rec_x_OMP_mod_v1 = x_OMP_mod_v1(:,1);
                end

                if ~isequal(rec_x_OMP_mod_v1, Path_arr)
                    error_count_OMP_mod_v1 = error_count_OMP_mod_v1 +1; % Increment count
                    flag=0;
                else
                    flag=1; % Correct result
                end 
                size_v1=size(x_OMP_mod_v1,2);
            end
%-------------------------------------------------------------------------%
%------------------ Recovery using PL-OMP & two missing link -----------------%
            if PLomp2
                ts_PLomp2 = tic;
                x_PLomp2 = PL_OMP2(n,h,y,Ar,L,B,C,temp_e,temp_f,temp_g,temp_dest,src.Node_id,dest.Node_id,'SE');
                te_PLomp2(pkt_i) = toc(ts_PLomp2) ;

                if isempty(x_PLomp2)
                    rec_x_PLomp2 = zeros(E,1);
                else
                    rec_x_PLomp2 = x_PLomp2(:,1);
                end

                if ~isequal(rec_x_PLomp2, Path_arr)
                    error_count_PLOMP2 = error_count_PLOMP2 +1; % Increment count
                    if flag==1
                        fprintf("BUG here\n")
                    end
                end 
                size_PLomp2=size(x_PLomp2,2);
            end
%-------------------------------------------------------------------------%
            
            
%------------------ Recovery using modified gOMP -------------------------%
            
            %N = 2;%floor(min(h,m(m_index)/h));

%-------------------- PL-gOMP with one missing link ----------------------%
            if PLgomp
                gOMP_ts_module1 = tic;
                [x_gOMP_mod_v1,bad_fraction(pkt_i)] = gOMP_list_module1(h,y,Ar,N,n,temp_e,temp_f,temp_g,temp_dest,B,C,src.Node_id);   
                x_gOMP_mod_v1(abs(x_gOMP_mod_v1)<=0.001)=0;
                x_gOMP_mod_v1(abs(x_gOMP_mod_v1)>0.001)=1;
                orig_x_gOMP_mod_v1 = x_gOMP_mod_v1;

                jj=[];              % Invalid path indices from set of solutions
                for j=1:size(x_gOMP_mod_v1,2)
                        x_nminus1_sq = x_gOMP_mod_v1(:,j);
                        x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
                        opt = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
                        if opt==0
                            jj=[jj j];
                        end
                end
                x_gOMP_mod_v1(:,jj)=[];
                gOMP_te_module1(pkt_i) = toc(gOMP_ts_module1) ;
                
                % Arrage in ascending residual norm
                residual_norm=zeros(1,size(x_gOMP_mod_v1,2));
                for j=1:size(x_gOMP_mod_v1,2)
                    final_residue = y - Ar*x_gOMP_mod_v1(:,j);
                    residual_norm(j) = norm(final_residue);
                end
                [~,order]= sort(residual_norm);
                x_gOMP_mod_v1 = x_gOMP_mod_v1(:,order);

                if isempty(x_gOMP_mod_v1)
                    rec_x_gOMP_mod_v1 = zeros(E,1);
                    flag_v1=0;
                else
                    rec_x_gOMP_mod_v1 = x_gOMP_mod_v1(:,1);
                    flag_v1=1;
                end
                if ~isequal(rec_x_gOMP_mod_v1, Path_arr)
                   error_count_gOMP_mod_v1 = error_count_gOMP_mod_v1 +1; % Increment count
                end
            end
%}
%------------------------------- L-gOMP  ---------------------------------%
            if Lgomp
                gOMP_ts_module2 = tic;
                x_gOMP_mod_v2 = gOMP_list_module2(h,y,Ar,N);         % List at first level only
                x_gOMP_mod_v2(abs(x_gOMP_mod_v2)<=0.001)=0;
                x_gOMP_mod_v2(abs(x_gOMP_mod_v2)>0.001)=1;
                orig_x_gOMP_mod_v2 = x_gOMP_mod_v2;

                jj=[];              % Invalid path indices from set of solutions
                for j=1:size(x_gOMP_mod_v2,2)
                        x_nminus1_sq = x_gOMP_mod_v2(:,j);
                        x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
                        opt = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
                        if opt==0
                            jj=[jj j];
                        end
                end
                x_gOMP_mod_v2(:,jj)=[];
                gOMP_te_module2(pkt_i) = toc(gOMP_ts_module2) ;
                
                % Arrange in ascending residual norm
                residual_norm=zeros(1,size(x_gOMP_mod_v2,2));
                for j=1:size(x_gOMP_mod_v2,2)
                    final_residue = y - Ar*x_gOMP_mod_v2(:,j);
                    residual_norm(j) = norm(final_residue);
                end
                [~,order]= sort(residual_norm);
                x_gOMP_mod_v2 = x_gOMP_mod_v2(:,order);

                if isempty(x_gOMP_mod_v2)
                    rec_x_gOMP_mod_v2 = zeros(E,1);
                    flag_v2=0;
                else
                    rec_x_gOMP_mod_v2 = x_gOMP_mod_v2(:,1);
                    flag_v2=1;
                end

                if ~isequal(rec_x_gOMP_mod_v2, Path_arr)
                   error_count_gOMP_mod_v2 = error_count_gOMP_mod_v2 +1; % Increment count
                   if ismember(Path_arr',x_gOMP_mod_v2','rows')
                        fprintf("RESIDUE BUG gOMP\n") 
                   end

                end
            end
            %{
            [flag_opt, flag_id]= ismember(rec_x_gOMP_mod_v2',x_gOMP_mod_v1','rows') ;
            if ~flag_opt
                if flag_v2==1   
                    fprintf("BUG") %only when module 2 itself has incorrect result
                    uiwait
                else
                    fprintf("ignore")
                end
            else
                if flag_v2==1 && flag_v1==0
                    fprintf("which column m2 not m1? %d\n",flag_id);
                    %uiwait
                elseif flag_v2==1 && flag_v1==1
                fprintf("which column? m2 m1 %d\n",flag_id);
                %uiwait
                end
            end
            
%-------------------------------------------------------------------------%            
            
%------------ Recovery using OMP when topology is known ------------------%
            % Selecting columns of Ar using topology adjacency matrix
            temp_k=[];
            kk=1;
            for i=1:(n-1)
                for j=1:n
                    if i==j
                    elseif graph_adj_mat(i,j)==1
                        temp_k=[temp_k kk];
                        kk=kk+1;
                    else
                        kk=kk+1;
                    end
                end
            end
            A_topo = Ar(:,temp_k);
            x_topo = OMP_modified(h,y,A_topo,L);
            x_topo(abs(x_topo)<=0.001)=0;
            x_topo(abs(x_topo)>0.001)=1;
            
            % Embed zeros for x_topo to x_(n-1)^2
            
            x_OMP_top=zeros(E,size(x_topo,2));
            x_check = x_OMP_top;
            temp_check = temp_k;
            for uu=1:size(x_topo,2)
                x_check(temp_check,uu) = x_topo(:,uu); 
                kk=1;ki=1;  % Count variables
                for i=1:(n-1)
                    for j=1:n
                        if i==j
                        elseif graph_adj_mat(i,j)==1
                            x_OMP_top(kk,uu)=x_topo(ki,uu);
                            kk=kk+1;
                            ki=ki+1;
                        else
                             x_OMP_top(kk,uu)=0;
                             kk=kk+1;
                        end
                    end
                end
            end
            if ~isequal(x_check, x_OMP_top)
                fprtintf("Bug")
                break
            end
            orig_x_OMP_top = x_OMP_top;
            ji=[];
            for j=1:size(x_OMP_top,2)
                x_nminus1_sq = x_OMP_top(:,j);
                x_n_sq = embed_extra_zeros(n,x_nminus1_sq,temp_e,temp_f,temp_g,temp_dest);
                opt1 = Adjacency_constraint(x_n_sq,B,C,h,src.Node_id);
              % opt2=valid_path(x_OMP_top(:,j),n,h,src.Node_id);
                if opt1==0
                     ji=[ji j];
                end    
            end
            x_OMP_top(:,ji)=[];
           
            if isempty(x_OMP_top)
                rec_x_OMP_top = zeros(E,1);%orig_x_OMP_top(:,1);
            else
                rec_x_OMP_top = x_OMP_top(:,1);
            end
            
            if ~isequal(rec_x_OMP_top, Path_arr)
               error_count_OMP_topo = error_count_OMP_topo +1; % Increment count
            end
           %} 
           %}
%------------------------------------ CVX --------------------------------%
            %x = cvx_solver(E,y,Ar);
            
            if cvx
                x = zeros(E,1);
                cvx_ts_module = tic;
                cvx_begin quiet
                cvx_precision default
                variable x(E) %binary  
                minimize( norm(x, 1 ) )
                subject to
                    Ar*x == y;

                cvx_end

                x(abs(x)<=0.001)=0;
                x(abs(x)>0.001)=1;
                cvx_te_module(pkt_i) = toc(cvx_ts_module) ;

                % Comparing recovered path wth original path
                if ~isequal(x, Path_arr)
                    error_count = error_count+1; %Increment if path recovered is different from path travelled
                end
            end
             %}
%-------------------------------------------------------------------------%
            % HALT CONDITION
            if error_count_OMP_mod_v2 == error_threshold % when threshold reached
                break
            end
        end
        
        % Time 
        if PLomp
            OMP_t_module1(m_index) = mean(OMP_te_module1);
        end
        if PLomp2
            t_PLomp2(m_index) = mean(te_PLomp2);
        end
        if Lomp
            OMP_t_module2(m_index) = mean(OMP_te_module2);
        end
        if PLgomp
            gOMP_t_module1(m_index) = mean(gOMP_te_module1);
        end
        if Lgomp
            gOMP_t_module2(m_index) = mean(gOMP_te_module2);
        end
        if cvx
            cvx_t_module(m_index) = mean(cvx_te_module);
        end
        %t_module3(m_index) = mean(te_module3);
        %t_Brute(m_index) = mean(te_brute);
        %fprintf("Time complexity:") ;
        %display(OMP_t_module1);
        
        
        % List size
        %fprintf("size PL-gomp list:");display(size(orig_x_gOMP_mod_v1,2));
        %fprintf("size L-gomp list:");display(size(orig_x_gOMP_mod_v2,2));
        %fprintf("size PL-omp list:");display(size_v1) ;
        %fprintf("size L-omp list:");display(size_v2);
        
        %}
        %bad_fraction_overall(m_index) = mean(bad_fraction); 
        %fprintf("Fraction of bad solutions in gOMP module 1:")
        %display(bad_fraction_overall);
        %list_wo_path_rate(m_index)= list_wo_path_count/pkt_count;
        
        error_rate(m_index) = error_count/pkt_count;
        error_rate_OMP(m_index) = error_count_OMP/pkt_count;
        error_rate_OMP_mod_v1(m_index) = error_count_OMP_mod_v1/pkt_count;
        error_rate_PLOMP2(m_index)=  error_count_PLOMP2/pkt_count;
        error_rate_OMP_mod_v2(m_index) =  error_count_OMP_mod_v2/pkt_count;
        %error_rate_OMP_topo(m_index) = error_count_OMP_topo/pkt_count;
        %error_rate_OMP_mod(m_index) =  error_count_OMP_mod/pkt_count;
        
        error_rate_gOMP(m_index) =  error_count_gOMP/pkt_count;
        error_rate_gOMP_mod_v1(m_index) =  error_count_gOMP_mod_v1/pkt_count;
        error_rate_gOMP_mod_v2(m_index)= error_count_gOMP_mod_v2/pkt_count;
        %error_rate_stOMP(m_index) = error_count_stOMP/pkt_count;
        %error_rate_cosamp(m_index) = error_count_cosamp/pkt_count;
        %}
        % DISPLAY Error rates
        %fprintf("Error Rate:%f for column size %d\n",error_rate(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size %d OMP\n",error_rate_OMP(m_index),m(m_index));
        fprintf("Error Rate:%f for column size %d PL-OMP\n",error_rate_OMP_mod_v1(m_index),m(m_index));
        fprintf("Error Rate:%f for column size %d PL-OMP 2\n",error_rate_PLOMP2(m_index),m(m_index));
        fprintf("Error Rate:%f for column size %d L-OMP\n",error_rate_OMP_mod_v2(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size %d mod OMP \n",error_rate_OMP_mod(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size %d OMP topo\n",error_rate_OMP_topo(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size %d list wo path\n",list_wo_path_rate(m_index),m(m_index));
      
        %fprintf("Error Rate:%f for column size %d gOMP \n",error_rate_gOMP(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size %d PL-gOMP\n",error_rate_gOMP_mod_v1(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size %d L-gOMP\n",error_rate_gOMP_mod_v2(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size %d stomp\n",error_rate_stOMP(m_index),m(m_index));
        %fprintf("Error Rate:%f for column size %d cosamp\n",error_rate_cosamp(m_index),m(m_index));
%}    
    end
    
      
    
    
%% PLOTS    
    
    %
    figure(2)
    semilogy(m,error_rate, 'yo-', 'LineWidth', 1.5);
    hold on
    semilogy(m,error_rate_OMP, 'b+-', 'LineWidth', 1.5);
    semilogy(m,error_rate_OMP_mod_v1, 'ms-', 'LineWidth', 1.5);
    semilogy(m,error_rate_OMP_mod_v2, 'gs-', 'LineWidth', 1.5);
    semilogy(m,error_rate_gOMP, 'k+-', 'LineWidth', 1.5);
    semilogy(m,error_rate_gOMP_mod_v1, 'cs-', 'LineWidth', 1.5);
    semilogy(m,error_rate_gOMP_mod_v2, 'rs-', 'LineWidth', 1.5);
    %semilogy(m,error_rate_stOMP, 'yo-', 'LineWidth', 2);
    %semilogy(m,error_rate_cosamp, 'ko-', 'LineWidth', 2);
    axis([0 100 0 1]);
    grid on
    legend('CVX','OMP','PL-OMP','L-OMP','gOMP','PL-gOMP','L-gOMP');
    title('Error rate vs column size SE');
    xlabel('Column size, m');
    ylabel('Error rate')
    hold off
    %}
%end

infos ={'Results for SEE:(I)CVX (II)1. OMP 2. PL-OMP 3. L-OMP (III)1. gOMP 2.PL-gOMP 3. L-gOMP';
   'n= no. of nodes';'h = path hop length/sparsity';'E= no. of possible edges';
   'm = column size of Ar for OMP';'graph_adj_mat = topology adjacency matrix'};
if save_result
   save SEE_n9_h8_fix.mat infos n E h graph_adj_mat path m error_rate error_rate_OMP error_rate_OMP_mod_v1 error_rate_OMP_mod_v2 error_rate_gOMP error_rate_gOMP_mod_v1 error_rate_gOMP_mod_v2 bad_fraction_overall OMP_t_module1 OMP_t_module2 gOMP_t_module1 gOMP_t_module2 cvx_t_module
end
