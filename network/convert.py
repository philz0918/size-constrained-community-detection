class convertData():

    def __init__(self, filename):

        self.input_file = filename

    def decreaseBy1(self, progress):

        file = open(self.input_file, "r")

        pair_dict = {}
        nodes = []
        comms = []

        for line in file:

            node = line.split()
            s_node = int(node[0]) - 1
            s_comm = int(node[1]) - 1

            if progress == "ntwk":

                nodes.append(s_node)
                comms.append(s_comm)
                pair_dict[s_node] = s_comm

            elif progress == "comm":

                if s_comm not in pair_dict.keys():

                    pair_dict[s_comm] = []
                    pair_dict[s_comm].append(s_node)

                else:
                    pair_dict[s_comm].append(s_node)

                nodes.append(s_node)
                comms.append(s_comm)

        return nodes, comms, pair_dict

    def genNewList(self, nodes, edges):
        n = 1
        newList = []
        for i in range(len(nodes)):
            if n < len(nodes):
                new = str(nodes[i]) + " " + str(edges[i]) + "\n"
                newList.append(new)
                n += 1
            else:
                new = str(nodes[i]) + " " + str(edges[i])
                newList.append(new)
                n += 1

        return newList

    def genCommNode(self, comm_vals):
        comm_list = []
        for comms in comm_vals:

            new = ""
            idx = 1
            for i in range(len(comms)):

                if idx < len(comms):
                    new += str(comms[i]) + " "
                    idx += 1
                else:
                    new += str(comms[i])

            # print(new)
            new += "\n"
            comm_list.append(new)

        return comm_list

    def savingNewG(self, root, txt_name, extension, newList):
        #need to change here
        newG = open(root + txt_name + "." + extension, "w")
        text_list = newList
        newG.writelines(text_list)

        return print("Processed!")



