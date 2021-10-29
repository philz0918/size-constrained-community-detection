
module helper

    using  PyCall
    using OrderedCollections
    ig = pyimport("igraph")

    mutable struct HelperStruct
        
        Q::Float64
        temp::Float64
        beta::Float64
        min::Int64
        max::Int64
        n::Int64
        
        G 
        comm::OrderedDict
        adj_list::Vector
    end

    #init variable struct
    function importing(g_file::String,g_min::Int64,g_max::Int64,g_beta::Float64)

        G = ig.Graph.Read(g_file)
        n = G.vcount()
        min = g_min
        max = g_max
        beta = g_beta
        Q = 0.0
        temp = 0.0

        #julia start with 1
        adj_list = Any[]
        comm = Dict()
        adj_temp = G.get_adjlist();
        for adj in adj_temp
            attemp = [node+1 for node in adj]
            push!(adj_list, attemp)
        end
        
        Help = HelperStruct(Q, temp, beta, min, max, n, G, comm, adj_list)

        return Help
    end

    function get_modularity(g, grouping::Vector)
        Q =ig.VertexClustering(g, grouping).q
        return Q
    end

    #update or not
    function helper_check_move(old_modularity::Float64, new_modularity::Float64, temp::Float64)
        
        if new_modularity >= old_modularity
            return true
        end

        thres = exp(-abs(old_modularity -new_modularity)/temp)
        if rand() < thres
            return true
        else 
            return false
        end
    end

    #reset value from 1
    function helper_reset_comm_values(Help::HelperStruct,updated_comm::OrderedDict)
        
        comm_values= Set([idx for idx in values(updated_comm)]);
        comm_values= [idx for  idx in comm_values];
        new_values =collect(1:length(comm_values));
        mapping = Dict()
        
        for c in new_values 
            temp_mapping = Dict(comm_values[c] => new_values[c])
            merge!(mapping, temp_mapping)
        end
        
        new_comm =OrderedDict()
        for (node, comms) in updated_comm
            temp_comm = Dict(node => mapping[comms])
            merge!(new_comm, temp_comm)
        end
        Help.comm = new_comm
        
    end

    #add beta regularization
    function helper_get_beta_modularity(Help::HelperStruct,b_comm::Vector)

        comms_grouped = [[node for (node, c) in b_comm if c ==i] for i in Set([val for val in values(b_comm)])];
        sizes = [size(x)[1] for x in  comms_grouped]
        reg_temp = 0 
        #println(sizes)
        for size in sizes 
            if size < min
                reg_temp += (min - size)^2;
            elseif size > max
                reg_temp += (size - max)^2;
            end
        end
        reg_param = beta*reg_temp;

        return reg_param
    end

    #check number of each groups
    function constraint(Help::HelperStruct,list::Vector)
        @eval import StatsBase: countmap
        cntmap = countmap(list);
        cnt_val = Set([val for val in values(cntmap)]);
        #println([val for val in values(cntmap)])
        for cnt in cnt_val 
            
            if  cnt < Help.min
                return false
            elseif cnt > Help.max 
                return false
            else 
                continue
            end
        end
        return true
    end   

end

module simulator

    using Random
    using ProgressBars
    using Combinatorics
    import ..helper
    import StatsBase: countmap
    using OrderedCollections

    g_modularity = 0

    #initial phase
    function init_graph(g_file::String,min::Int64,max::Int64,beta::Float64)
        @eval using OrderedCollections
        
        Help = helper.importing(g_file,min,max,beta);
        #init_num_groups = Int64(floor(nodes/min));
        init_num_groups = Int64(floor(Help.n/Help.min));
        perm = randperm(Help.n);

        init_partition = Int64[];
        for i in 1:Help.n
            val = perm[i]%init_num_groups;
            append!(init_partition,val);
        end
        
        comm = OrderedDict(i => init_partition[i] for i in 1:Help.n);
        comm_list = [group for group in values(comm)];
        
        modularity = helper.get_modularity(Help.G, comm_list);

        #assign global variables; modularity, comm

        Help.comm = comm 
        Help.Q = modularity


        return Help
    end

    function local_move(Help::helper.HelperStruct)

        #inherit HelperStruct from helper
        l_comm = deepcopy(Help.comm);
        adj = Help.adj_list;

        #for debugging log
        tup_list =  Any[];
        append!(tup_list,[(Help.Q,"      ",[val for val in values(l_comm)])]);

        ## for while loop condition
        global num_check = false;

        while num_check == false
            adj_check = true;
            #=
            this while loop for checking the community that trying to move contains connected node
            =#
            while adj_check == true
                global selected_comm_idx = rand(Set([idx for idx in values(l_comm)]));
                global selected_node = rand([node for (node, comm_idx) in l_comm if comm_idx == selected_comm_idx]);

                #adjacent nodes of selected nodes
                adj_nodes = adj[selected_node];
                global adj_selected_comm_idx = Set([l_comm[n] for n in adj_nodes]);
                
                if length(adj_selected_comm_idx) == 1 && selected_comm_idx in adj_selected_comm_idx
                    
                    adj_check = true;
                    #adj_check = false;
                else 
                    adj_check = false;
                end
            end
        
            #println("End here")
            listy = [val for val in adj_selected_comm_idx if  val != selected_comm_idx]
            #=
            This is the only case when max == n , we can extnd new community
            =#
            if length(Set([idx for idx in values(l_comm)])) == 1
                append!(listy, maximum(values(l_comm))+1);
            end

            global new_selected_comm_idx = rand(Set(listy));

            #update temporary community 
            temp_comm = copy(l_comm)
        
            global temp_comm[selected_node] = new_selected_comm_idx
            global updated_comm_list =[val for val in values(temp_comm)]
        
            #check constraint of minmax
            #=
            if Help.temp > 0.00001
                num_check = true
            else 
                num_check = helper.constraint(Help,updated_comm_list)
                if num_check == false
                    println(num_check)
                end
            end
            =#
            
            num_check = helper.constraint(Help,updated_comm_list)
            
            if num_check == true
                break
            end
            
        end
        #debugging log
        append!(tup_list, [(num_check, "change ", selected_node, selected_comm_idx, new_selected_comm_idx, updated_comm_list)]);

        modularity = helper.get_modularity(Help.G, updated_comm_list);
        update_check = helper.helper_check_move(Help.Q, modularity, Help.temp);

        if update_check == true

            Help.Q = modularity;
            temp_comm  = sort!(OrderedDict(temp_comm), byvalue = false);
            Help.comm = temp_comm;

        end
        
        return tup_list
    end

    function merge(Help::helper.HelperStruct)
        @eval using Combinatorics
        community= copy(Help.comm)
        set_comm = Set([group for group in values(community)])
        #check the case when it has only 1 community
        if length(set_comm) == 1
            return community
        end
        
        group_count = OrderedDict(countmap([group for group in values(community)]))
        combi = collect(combinations(minimum(set_comm):maximum(set_comm),2))
    
        combi_all = Any[]
        for com in combi
        
            sum = group_count[com[1]]+group_count[com[2]]
            if Help.min <= sum <= Help.max 
                push!(combi_all, com)
            end
        end
        
        #if there is no satisfied combination 
        if length(combi_all) == 0
            return community
        end

        # choose combination containing (source, target)
        selected_combi = rand(combi_all)
        source_nodes = [nodes for (nodes, comm) in community if comm == selected_combi[1]]
        # source -> target comms
        for n in source_nodes
            community[n] = selected_combi[2]
        end
        

        return community
    end

    function split(Help::helper.HelperStruct)
        @eval import StatsBase: countmap
        community = copy(Help.comm)

        group_count = countmap([group for group in values(community)])
        cnt_comm =[val for val in Set(values(group_count)) if val >= Help.min*2]
        
        #check the condition for elligible for split group has at least min*2 
        if length(cnt_comm) ==0
            return community
        end
        
        candidate_group = [group for (group, num) in group_count if num in cnt_comm]
        selected_comm = rand(candidate_group)
        selected_nodes = [node for (node,comm) in community if comm == selected_comm]
        slicenum = rand(collect(Help.min:group_count[selected_comm]-Help.min))
        
        #=
        split groups, length(group_count) => total number of groups
        new_group_nodes => nodes being changed
        =#

        new_group_nodes = selected_nodes[1:slicenum]
        new_comm = length(group_count)+1

        for n in new_group_nodes
            community[n] = new_comm
        end

        return community
    end


    function globalmove(Help::helper.HelperStruct, p_m::Float64 = 0.495, p_s::Float64 = 0.495)

        p = rand()
        
        if p <= p_m
            community = merge(Help)
            move = "M"
            #println(move)
        elseif p_m <= p <= p_m+p_s
            community = split(Help)
            move = "S"
            #println(move)
        else
            local_move(Help)
            move = "L"
            #println(move)
        end
        
        #different for L
        #updated_comm_list = [val for val in values(community)]

        if move == "M" || move == "S"
            updated_comm_list = [val for val in values(community)]
            modularity = helper.get_modularity(Help.G, updated_comm_list);
            update_check = helper.helper_check_move(Help.Q, modularity, Help.temp);

            if update_check == true
    
                Help.Q = modularity;
                temp_comm  = sort!(OrderedDict(community), byvalue = false);
                Help.comm = temp_comm;

            end
        end
        return 
    end

    function annealing(g_file::String,min::Int64,max::Int64,temp::Float64,stop_epsilon::Float64, cooling_factor::Float64,beta::Float64,name::String)
        @eval using ProgressBars
        @eval using DelimitedFiles
        Help = init_graph(g_file,min,max,beta);
        Help.temp = temp
        
        num_moves = floor(Int,(log(stop_epsilon/temp)/log(cooling_factor))+1);
        total_list =  Any[]
        q_list = Any[]
        
        for i in (1:num_moves)
            #println(helper.temp)

            for _ in (1:Help.n)
                globalmove(Help);
                helper.helper_reset_comm_values(Help, Help.comm);
            end
                                    
            for _ in (1:Help.n^2)
                t_list = local_move(Help);
                helper.helper_reset_comm_values(Help,Help.comm);
                append!(total_list, t_list)
            end
            
            num_check = helper.constraint(Help,updated_comm_list)
                
            #println("     ")
            println(Help.Q)
            append!(q_list, Help.Q)
            temp *= cooling_factor;
            Help.temp = temp
        end
        
        #=
        open("result/tuple_list_"* name *".txt", "w") do file
            writedlm(file, total_list)
        end
        =#
        open("result/q_list_"* name *".txt", "w") do file
            writedlm(file, q_list)
        end
        
        return (Help.comm, Help.Q)
    end

end