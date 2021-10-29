import igraph as ig
from convert import convertData

"""
benchmark : if given network is benchmark network which require to be processed
root : saving directory
ntwk_root : benchmark network file
ntwk : new network file name
ntwk_fm : saving network format default(igraph) : "dat", networkX : "edgefile"
coom_root : benchmark community file 
comm : new community file name
network : network file(.gml) that does not require further processing 
"""


def Network (benchmark = False ,root = None , ntwk_root = None, ntwk = None, ntwk_fm = "dat",comm_root = None, comm = None, network = None, gt_comm = None) :

    if benchmark is True :

        if ntwk_root is None :
            pass
        else :
            '''
            for network from benchmark
            input : (.dat) pairs file
            output : saving with pairs (decreased by 1)
            '''

            convert = convertData(ntwk_root)
            nodes, comms,_ = convert.decreaseBy1("ntwk")
            convert_list = convert.genNewList(nodes, comms)
            convert.savingNewG(root,ntwk,ntwk_fm,convert_list)

        if comm_root is None:
            pass
        else :
            '''
            for community from benchmark 
            input : community.dat
            output : list of community and community files that community aggregated by nodes (Using for NMI)
              (ex) 0 1 2 4 5 20 39 : these nodes are in same community 
            '''

            convert_comm = convertData(comm_root)
            _, comm_num, pair_dict = convert_comm.decreaseBy1("comm")
            comm_list = convert_comm.genCommNode(list(pair_dict.values()))
            convert_comm.savingNewG(root,comm,"dat",comm_list)

        preG = ig.read(root+ntwk+".dat")
        preG = preG.as_undirected()
        preG.save(root+ntwk+".gml", format = "gml")
        print("the number of communities in this network is {}".format(len(set(comm_num))))
        true_comm = len(set(comm_num))

        # get ground truth modularity of benchmark network
        G = ig.read(root+ntwk+".gml")
        partition = ig.VertexClustering(G, membership=comm_num)
        optim_modul = partition.q

        #ground truth
        print("ground truth modularity {}".format(optim_modul))



    else :

        try :
            assert network is None, "Provide network file(.gml)"
        except AssertionError as msg :
            print(msg)

        G = ig.read(network)
        print("the number of communities in this network is {}".format(gt_comm))
        true_comm = gt_comm

    return G, true_comm


if __name__ == "__main__" :
    '''
    for n in [500,1000] :
        for mu in [0.1,0.3,0.75] :
            print("{} nodes, {} mu".format(n, mu))
            G, true_comm = Network(benchmark= True,
                        root = "/Users/macbookpro/Constrained_modularity/python/generated/",
                        ntwk_root= "/Users/macbookpro/Constrained_modularity/python/benchmark/network_"+str(n)+"_"+str(mu)+".dat",
                        ntwk = "ntwk"+str(n)+"_"+str(mu),
                        comm_root = "/Users/macbookpro/Constrained_modularity/python/benchmark/community_"+str(n)+"_"+str(mu)+".dat",
                        comm = "comm"+str(n)+"_"+str(mu))
    '''
    G, true_comm = Network(benchmark=True,
                           root="/Users/macbookpro/Constrained_modularity/python/generated/",
                           ntwk_root="/Users/macbookpro/Constrained_modularity/python/ntwkfile/network_50k_0.05-mu_1.dat",
                           ntwk="ntwk_test_1",
                           comm_root="/Users/macbookpro/Constrained_modularity/python/ntwkfile/community_50k_0.05-mu_1.dat",
                           comm="comm_test_1")
