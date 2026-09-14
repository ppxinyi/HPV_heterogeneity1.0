from tree_class import *
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import collections  as mc
from collections import Counter
import pylab


#####open pickle junctionRes
#junctionListDic: diretory of junctionList
def open_junctionList(junctionListDicPath):
    with open(junctionListDicPath,'rb') as inputFile:
        junctionListDic = pickle.load(inputFile)
    return junctionListDic

######generate dictionary of junction local depth
#extract depth Dictionary for each junction
# def depth_dic(junctionListDic):
#     depthDic = {}
#     for read in junctionListDic:
#         for info in junctionListDic[read]:
#             i = 0
#             for each in info['chro']:
#                 parent_node = each + ':' + str(info['junction'][i]) 
#                 junction = parent_node
#                 depth = info['depth'][i]
#                 if junction in depthDic:
#                     depthDic[junction].append(depth)
#                 else:
#                     depthDic[junction] = [depth]
#                 i += 1
#     newDepthDic = {}
#     for junction in depthDic:
#         newDepthDic[junction] = max(depthDic[junction])
#     return newDepthDic

######generate dictionary of junction
#readNameList: readNameList start to search
def abundent_dic(readNameList,junctionListDic):
    abundentDic = {}
    tempDic = {}
    for read in readNameList:
        for info in junctionListDic[read]:
            # if read == '436e9433-2402-4f77-9421-f203875eb89f':
            i = 0
            for each in info['chro']:
                
                parent_node = each.split(':')[0] + ':' + info['junction'][i].split(':')[0]
                node = each.split(':')[1] + ':' + info['junction'][i].split(':')[1]
                if read in abundentDic:
                    abundentDic[read].append(parent_node + ';' + node + ';' + str(i)) #store the order of the node
                else:
                    abundentDic[read] = [parent_node + ';' + node + ';' + str(i)]
                i += 1
                    
    return abundentDic

######generate dictionary of junction, key = readName
#abundentDic: return from abundent_dic
def junction_dic(abundentDic):
    tempDic= {}
    for read in abundentDic:
        for each in abundentDic[read]:
            junction = ';'.join(each.split(';')[:-1])
            if junction in tempDic:
                tempDic[junction].append(read)
            else:
                tempDic[junction] = [read]
    maxNum = 0
    print(abundentDic)
    print(tempDic)            
    for junction in tempDic:
        if len(tempDic[junction]) > maxNum:
            maxNum = len(tempDic[junction])
    newReadNameList = []
    print(maxNum)
    print(tempDic)
    print(abundentDic)
    targetJunction = None
    for junction in tempDic:
        if len(tempDic[junction]) == maxNum:
            newReadNameList = tempDic[junction]
            targetJunction = junction
            break
    newReadNameList = list(set(newReadNameList))
    return tempDic,newReadNameList,targetJunction

######search more most abundent junction
#readNameList: readNameList start to search
def search_most_junction(readNameList,junctionListDic):
    abundentDic = abundent_dic(readNameList,junctionListDic)
    ##### debug chr2####
    target = "2:151136063;gi|333031|lcl|HPV16REF.1|:6590"
    owners = [r for r, nodes in abundentDic.items() if any(target in x for x in nodes)]
    print("[DBG] target owners in step①:", owners)
    ####3 degure chr2######### 
    tempDic,newReadNameList,targetJunction = junction_dic(abundentDic)            
    newAbundentDic = abundent_dic(newReadNameList,junctionListDic)
    #####chr2 ####
    print("[DBG] target present after step③?",
          any(any(target in x for x in nodes) for nodes in newAbundentDic.values()))
    ##### chr2 ########
    return newAbundentDic,targetJunction,newReadNameList


######search target junction in read
def search_junction(readName,targetJunction,abundentDic):
    position = 0
    i = 0
    has = False
    for each in abundentDic[readName]:
        node = ';'.join(each.split(';')[:-1])
        if i < len(abundentDic[readName]):
            if node == targetJunction:
#                 print(readName,i)
                position = i
                has = True
                #print(position)
                break
        i += 1
    return has,i

######search target junction in read
def search_previous_junction(readName,targetJunction,abundentDic):
    preNode = None
    i = 0
    for each in abundentDic[readName]:
        node = ';'.join(each.split(';')[:-1])
        if i > 0:
            if node == targetJunction:
                preNode = abundentDic[readName][i-1]
        i += 1
    return preNode

######store junction for each read
def store_next_junction_each_read(tree,readName,abundentDic,targetJunction):
    rootNode = targetJunction
    has,position = search_junction(readName,targetJunction,abundentDic)
    if has:
        i = position
        for each in abundentDic[readName][position:]:
            node = ';'.join(each.split(';')[:-1])
            order_of_node = each.split(';')[-1]
            #print(node)
            if i < len(abundentDic[readName])-1:
                #check whether nextNode exist in previous tree, if so, it's duplication
                next_node = ';'.join(abundentDic[readName][i+1].split(';')[:-1])
                read_name_in_node = readName + ';' + order_of_node
                read_name_in_next_node = readName + ';' + abundentDic[readName][i+1].split(';')[-1]
                if tree.breadth_first_search(rootNode,node) == True and tree.breadth_first_search(rootNode,next_node) == True:
                    # if abundentDic[readName][i+1] == 'gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:1':
                    #     print('\n')
                    #     print('node',node)
                    #     print('\n')
                    tree.add_readName_in_node(next_node,read_name_in_next_node,node,read_name_in_node)
                    tree.add_endNode_in_node(next_node,read_name_in_next_node,node,read_name_in_node)
                else:
                    tree.add_node(next_node,read_name_in_next_node,node,read_name_in_node)
                    # if abundentDic[readName][i+1] == 'gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:1':
                    #     print('\n')
                    #     print('node',node)
                    #     print('\n')
            
            i += 1
    return tree

######store junction for each read
def store_previous_junction_each_read(tree,readName,abundentDic,targetJunction):
    rootNode = targetJunction
    
    has,position = search_junction(readName,targetJunction,abundentDic)
    if has:
        i = position
        
        for i in range(position,0,-1):
            node = ';'.join(abundentDic[readName][i].split(';')[:-1])
            order_of_node = abundentDic[readName][i].split(';')[-1]

            #print(abundentDic[readName][i])
            #print(i,targetJunction,abundentDic[readName])
            #check whether preNode exist in previous tree, if so, it's duplication
            #XInyi change 1 to 0
            if i > 0:
                pre_node = ';'.join(abundentDic[readName][i-1].split(';')[:-1])
                read_name_in_node = readName + ';' + order_of_node
                read_name_in_pre_node = readName + ';' + abundentDic[readName][i-1].split(';')[-1]
                if tree.breadth_first_search(rootNode,node) == True and tree.breadth_first_search(rootNode,pre_node) == True:
                    #duplication
                    tree.add_readName_in_node(pre_node,read_name_in_pre_node,node,read_name_in_node)
                    tree.add_endNode_in_node(pre_node,read_name_in_pre_node,node,read_name_in_node)
                else:
                    tree.add_node(pre_node,read_name_in_pre_node,node,read_name_in_node)
            
            
    return tree

#####strore junction starting from this targetJunction
def store_forward_junction(readNameList,targetJunction,abundentDic):
    nextTree = Tree()
    #initial tree
    has,position = search_junction(readNameList[0],targetJunction,abundentDic)
    if has:
        order_of_node = abundentDic[readNameList[0]][position].split(';')[-1]
        read_name_in_node = readNameList[0] + ';' + order_of_node
        nextTree.add_node(targetJunction,read_name_in_node)

    for read in readNameList[1:]:
        has,position = search_junction(read,targetJunction,abundentDic)
        if has:
            order_of_node = abundentDic[read][position].split(';')[-1]
            read_name_in_node = read + ';' + order_of_node
            nextTree.add_readName_in_node(targetJunction,read_name_in_node)
    for read in readNameList:
        #print(read)
        store_next_junction_each_read(nextTree,read,abundentDic,targetJunction)
    return nextTree    

#####strore junction starting from this targetJunction
def store_backward_junction(readNameList,targetJunction,abundentDic):
    #initial tree
    nextTree = Tree()
    has,position = search_junction(readNameList[0],targetJunction,abundentDic)
    if has:
        order_of_node = abundentDic[readNameList[0]][position].split(';')[-1]
        read_name_in_node = readNameList[0] + ';' + order_of_node
        nextTree.add_node(targetJunction,read_name_in_node)
    for read in readNameList[1:]:
        has,position = search_junction(read,targetJunction,abundentDic)
        if has:
            order_of_node = abundentDic[read][position].split(';')[-1]
            read_name_in_node = read + ';' + order_of_node
            nextTree.add_readName_in_node(targetJunction,read_name_in_node)
    for read in readNameList:
        #print(read)
        store_previous_junction_each_read(nextTree,read,abundentDic,targetJunction)
        #nextTree.print_tree_breadth_first(targetJunction)
    return nextTree   

######iteratively assemble the tree
#junctionListDic
def build_tree(junctionListDic):
    junctionList = junctionListDic.keys()
    allList = junctionListDic.keys()
    allAssemblies = []
    while junctionList != []:
        abundentDic,targetJunction,readNameList = search_most_junction(junctionList,junctionListDic)
        if readNameList != [] and targetJunction != '' and abundentDic != {}:
            tree_forward = store_forward_junction(readNameList,targetJunction,abundentDic)
            resList_forward = []
        #     eachCluster = []
            # print("t:",targetJunction)
            tree_forward.walk_assembly(targetJunction,[],[],[],[],resList_forward)
        #     for read in resList_forward:
                #print(read)
        #     tree_backward = store_backward_junction(readNameList,targetJunction,abundentDic)
        #     resList_backward = []
        #     tree_backward.walk_assembly(targetJunction,[],[],[],[],resList_backward)
        #     eachCluster.append([read,resList_backward])
        #     if targetJunction in allAssemblies:        
        #         allAssemblies[targetJunction].append(eachCluster)
        #     else:
        #         allAssemblies[targetJunction] = [eachCluster]
            tree_backward = store_backward_junction(readNameList,targetJunction,abundentDic)
            resList_backward = []
            tree_backward.walk_assembly(targetJunction,[],[],[],[],resList_backward)
            
            allAssemblies.append([resList_forward,resList_backward])
            newJunctionList = []
        # print(len(set(abundentDic.keys())))
            for each in set(junctionList):
                if each not in set(abundentDic.keys()):
                    newJunctionList.append(each)
            junctionList = newJunctionList
        else:
            break

    return allAssemblies


#######filter low coverage paths
#allAssemblies: result returned from algorithms;
#cov: read coverage for this sample
#depthDic: local read depth dictionary
# def filter_read(allAssemblies,cov,depthDic):
#     filteredAllAssembiles = []
#     for direction in allAssemblies:
#         newdirection = []
#         for cluster in direction:
#             newCluster = []
#             for path in cluster:
                
#     #             filtered = False
#     #             for junction in path[1]:
#     #                 if len(path[1]) < 4 and junction[1] <= 4:
#     #                     filtered = True
#     #             i = len(path[1])-1
#     #             #print(len(path[1]))
#     #             for junction in path[1][::-1]:
#     #                 if junction[1] <= 4:
#     #                     i = i - 1
#                 i = len(path[0])-1
#                 #print(len(path[1]))
#                 for junction in path[0][::-1]:
#                     junction1 = junction.split(';')[0]
#                     junction2 = junction.split(';')[1]
#                     depth1 = depthDic[junction1]
#                     depth2 = depthDic[junction2]
#                     if min(depth1,depth2)<= cov/2:
#                         i = i-1
#                 #print(i)
#     #             if not filtered and not oneRead:
#                 newCluster.append([path[0][:i+1],path[1][:i+1],path[2][:i+1],path[3][:i+1]])
            
            
                

           
#             newdirection.append(newCluster)
#         filteredAllAssembiles.append(newdirection)

#     #remove all redundent path
#     uniqFilteredAllAssembiles = []
#     for direction in filteredAllAssembiles:
#         newdirection = []
#         for cluster in direction:
#             newCluster = []
#             allPath = []
#             for path in cluster:
#                 allPath.append(tuple(path[0]))
#             allPath = list(set(allPath))
#             for path1 in allPath:
#                 if path1 != tuple():
#                     for path2 in cluster:
#                         if path1 == tuple(path2[0]):
#                             newCluster.append(path2)
#                             break
#             #add forward and backward even if empty
#             newdirection.append(newCluster)
        
#         if direction != []:
#             uniqFilteredAllAssembiles.append(newdirection)
                    
#     return uniqFilteredAllAssembiles
                    
            

#####write the results to file
def write_assembly(allAssembiles,outputFile):
    with open(outputFile,'wb') as output:
        pickle.dump(allAssembiles,output)

        