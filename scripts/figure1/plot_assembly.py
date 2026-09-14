from tree_class import *
from store_junction import *
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import collections  as mc
from collections import Counter
import pylab

v_color = (194/235,207/235,162/235)
h_color = (164/235,121/235,158/235)

cm = pylab.get_cmap('tab20')

#######plot forward assemblies
#Parameters:
#uniqFilteredAllAssembiles: dictionary reported from filter_read()
#outputPath: path to save the plot
#sampleName: name of sample, used as prefix of the file name 
#Output:
#None
def plot_forward(uniqFilteredAllAssembiles,outputPath,sampleName):
    m = 0
    for path in uniqFilteredAllAssembiles:

        len_candidate = []
        lim_list = []
        #path[0] is forward tree
        #check whether forward is empty
        if path[0] != []:
            for candidate in path[0]:
                len_candidate.append(len(candidate[0]))
                lim = max(len(path[0]),max(len_candidate))
                lim_list.append(lim)
            
            

            nodeDic = {}
            nodeList = []
            for candidate in path[0]:
                for node in candidate[0]:
                    nodeList.append(node)
            if nodeList != []:
                fig,ax = plt.subplots(figsize = (lim*2,lim*2))
                n = 0
                for node in set(nodeList):
                    nodeDic[node] = cm(1.*n/(len(set(nodeList))))
                    # if node == 'gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:1':
                    #     print(n,1.*n/(len(set(nodeList))),cm(1.*n/(len(set(nodeList)))))
                    # if node == 'gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:3301':
                        # print(n,1.*n/(len(set(nodeList))),cm(1.*n/(len(set(nodeList)))))
                    n += 1
                # if 'gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:1' in nodeDic and 'gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:3301' in nodeDic:
                #     print(nodeDic['gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:1'],nodeDic['gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:3301'])


                j = 0
                for candidate in path[0]:
                    i = 0
                    for node in candidate[0]:
                        x = i+1
                        y = j+1
                        #add node
                        junction1_chro = node.split(';')[0].split(':')[0]
                        junction1_site = node.split(';')[0].split(':')[1]
                        if 'gi|333031|lcl|HPV16REF.1|:1' in node and 'gi|333031|lcl|HPV16REF.1|:3314' in node:
                            print(candidate[1][i],candidate[2][i])
                        if junction1_chro == 'gi|333031|lcl|HPV16REF.1|':
                            ax.add_patch(plt.Circle((x,y), radius=0.1,clip_on=False,color = nodeDic[node]))
                            site = node.split(';')[0].split(':')[1]
                            text = f'HPV:{site}'
                            label = 'V'
                        else:
                            ax.add_patch(plt.Rectangle((x-0.1,y-0.1), width = 0.2, height = 0.2,clip_on=False,color = nodeDic[node]))
                            site = node.split(';')[0].split(':')[1]
                            chro =node.split(';')[0].split(':')[0]
                            text = f'{chro}:{site}'
                            label = 'H'
                        ax.text(x, y+0.2, text, zorder=100,
                                ha='center', va='center', weight='bold',
                                style='italic',size = 8)

                        #add copyNumber
                        print(node,path[0][j][1][i])
                        copyNumber = 'X' + str(path[0][j][1][i][1])


                        ax.text(x, y+0.5, copyNumber, zorder=100,
                                ha='center', va='center', weight='bold',
                                style='italic',size = 8)
                        ax.text(x, y, label, zorder=100,
                                ha='center', va='center', weight='bold',
                                style='italic',size = 10)
                        offset = 0.2
                        junction2_chro = node.split(';')[1].split(':')[0]
                        junction2_site = node.split(';')[1].split(':')[1]
                        if junction2_chro == 'gi|333031|lcl|HPV16REF.1|':
                            x = i+1
                            y = j+1
                            ax.add_patch(plt.Circle((x+offset,y), radius=0.1,clip_on=False,color = nodeDic[node]))
                            site = node.split(';')[1].split(':')[1]
                            text = f'HPV:{site}'
                            label = 'V'
                        else:
                            ax.add_patch(plt.Rectangle((x-0.1+offset,y-0.1), width = 0.2, height = 0.2,clip_on=False,color = nodeDic[node]))
                            site = node.split(';')[1].split(':')[1]
                            chro =node.split(';')[1].split(':')[0]
                            text = f'{chro}:{site}'
                            label = 'H'
                        ax.text(x+offset, y-0.2, text, zorder=100,
                                ha='center', va='center', weight='bold',
                                style='italic',size = 8)
                        ax.text(x+offset, y, label, zorder=100,
                                ha='center', va='center', weight='bold',
                                style='italic',size = 10)
                        #add looping node
            #             loopingNode = str(path[0][j][3][i])
            #             if loopingNode != []:
            #                 print(node,loopingNode)
                        #add segment

                        if i < len(candidate[0])-1:
                            lines = []
                            junction1_chro_next = candidate[0][i+1].split(';')[0].split(':')[0]
                            #print(junction2_chro,junction1_chro_next)
                            if junction2_chro == 'gi|333031|lcl|HPV16REF.1|' and junction1_chro_next == 'gi|333031|lcl|HPV16REF.1|':

                                ax.plot([x+offset,x+1],[y,y],'-')
                                lines.append([(x+offset, y),(x+1,y)])
                            lc = mc.LineCollection(lines, linewidths=2,color = v_color)
                            ax.add_collection(lc)
                            lines = []
                            if junction2_chro != 'gi|333031|lcl|HPV16REF.1|' and junction1_chro_next != 'gi|333031|lcl|HPV16REF.1|':
                                ax.plot([x+offset,x+1],[y,y],'-')
                                lines.append([(x+offset, y),(x+1,y)])
                            lc = mc.LineCollection(lines, linewidths=2,color = h_color)
                            ax.add_collection(lc)

                        #check looping node
                        copyNumberLoop = []
                        for eachLoop in candidate[3][i]:
                            k = 0
                            #print(eachLoop)
                            for eachNode in candidate[0][:i]:
                                if ';'.join(eachLoop.split(';')[0:2])  == eachNode:
                                    #if 'gi|333031|lcl|HPV16REF.1|:1' in eachNode and 'gi|333031|lcl|HPV16REF.1|:3314' in eachNode:
                                    print(eachNode,i,k)
                                    print(candidate[0][i],candidate[0][k])
                                    #,eachLoop,eachLoop.split(';')[1],candidate[0],candidate[1])
                                    copyNumberLoop.append(str(i) + ';' + str(k))
                                k += 1
                        #plot the loops
                        countLoop = Counter(copyNumberLoop)
                        #print(countLoop)
                        lines = []
                        for eachLoop in countLoop:
                            offset = 1.1
                            start = int(eachLoop.split(';')[0])
                            end = int(eachLoop.split(';')[1])
                            copyNumber = countLoop[eachLoop]
                            #color
                            if candidate[0][end].split(';')[1] == 'gi|333031|lcl|HPV16REF.1|' and candidate[0][start].split(';')[0] == 'gi|333031|lcl|HPV16REF.1|':
                                loopColor = v_color
                            else:
                                loopColor = h_color

                            lines.append([(start+offset,y),(start+offset,y+0.3)])
                            lines.append([(start+offset,y+0.3),(end+offset,y+0.3)])
                            lines.append([(end+offset,y),(end+offset,y+0.3)])
                            ax.text((start+offset+end+offset)/2, y+0.4, 'X'+str(copyNumber), zorder=100,
                                ha='center', va='center', weight='bold',
                                style='italic',size = 8)
                        if lines != []:
                            lc = mc.LineCollection(lines, linewidths=2,color = loopColor)
                            ax.add_collection(lc)


                    #print(Counter(copyNumberLoop))




                        i += 1
                    j += 1


            ax.set_xlim(0,lim + 1)
            ax.set_ylim(0,lim + 1)
            fig.savefig(f'{outputPath}/{sampleName}.forward.{m}.png')
            plt.close(fig)
        m += 1

#######plot backward assemblies
#Parameters:
#uniqFilteredAllAssembiles: dictionary reported from filter_read()
#outputPath: path to save the plot
#sampleName: name of sample, used as prefix of the file name 
#Output:
#None
def plot_backward(uniqFilteredAllAssembiles,outputPath,sampleName):
    m = 0
    for path in uniqFilteredAllAssembiles:
        
        len_candidate = []
        lim_list = []
        #path[1] is back tree
        #check whether back is empty
        if path[1] != []:
            for candidate in path[1]:

                len_candidate.append(len(candidate[0]))
                lim = max(len(path[1]),max(len_candidate))
                lim_list.append(lim)
            fig,ax = plt.subplots(figsize = (lim*2,lim*2))

            nodeDic = {}
            nodeList = []
            for candidate in path[1]:
                for node in candidate[0]:
                    nodeList.append(node)
            n = 0
            for node in set(nodeList):
                nodeDic[node] = cm(1.*n/(len(set(nodeList))))
                # if '189883773' in node:
                #     print(node)
                #     print(n,1.*n/(len(set(nodeList))),cm(1.*n/(len(set(nodeList)))))
                n += 1
            j = 0
            for candidate in path[1]:
                i = 0
                for node in candidate[0]:
                    x = lim -i
                    y = j+1
                    #add node
                    junction1_chro = node.split(';')[1].split(':')[0]
                    junction1_site = node.split(';')[1].split(':')[1]
                    # if 'gi|333031|lcl|HPV16REF.1|:3720' in node and 'gi|333031|lcl|HPV16REF.1|:3314' in node:
                    #     print(candidate[1][i],candidate[2][i])
                    if junction1_chro == 'gi|333031|lcl|HPV16REF.1|':
                        ax.add_patch(plt.Circle((x,y), radius=0.1,clip_on=False,color = nodeDic[node]))
                        site = node.split(';')[0].split(':')[1]
                        text = f'HPV:{junction1_site}'
                        label = 'V'
                    else:
                        ax.add_patch(plt.Rectangle((x-0.1,y-0.1), width = 0.2, height = 0.2,clip_on=False,color = nodeDic[node]))
                        
                        text = f'{junction1_chro}:{junction1_site}'
                        
                        label = 'H'
                    ax.text(x, y+0.2, text, zorder=100,
                            ha='center', va='center', weight='bold',
                            style='italic',size = 8)

                    #add copyNumber
                    copyNumber = 'X' + str(path[1][j][1][i][1])


                    ax.text(x, y+0.5, copyNumber, zorder=100,
                            ha='center', va='center', weight='bold',
                            style='italic',size = 8)
                    ax.text(x, y, label, zorder=100,
                            ha='center', va='center', weight='bold',
                            style='italic',size = 10)
                    offset = 0.2
                    junction2_chro = node.split(';')[0].split(':')[0]
                    junction2_site = node.split(';')[0].split(':')[1]
                    if junction2_chro == 'gi|333031|lcl|HPV16REF.1|':
                        x = lim - i
                        y = j+1
                        ax.add_patch(plt.Circle((x-offset,y), radius=0.1,clip_on=False,color = nodeDic[node]))
                        text = f'HPV:{junction2_site}'
                        label = 'V'
                    else:
                        ax.add_patch(plt.Rectangle((x-0.1-offset,y-0.1), width = 0.2, height = 0.2,clip_on=False,color = nodeDic[node]))
                        text = f'{junction2_chro}:{junction2_site}'
                        label = 'H'
                    ax.text(x-offset, y-0.2, text, zorder=100,
                            ha='center', va='center', weight='bold',
                            style='italic',size = 8)
                    ax.text(x-offset, y, label, zorder=100,
                            ha='center', va='center', weight='bold',
                            style='italic',size = 10)
                    #add looping node
        #             loopingNode = str(path[0][j][3][i])
        #             if loopingNode != []:
        #                 print(node,loopingNode)
                    #add segment

                    if i < len(candidate[0])-1:
                        lines = []
                        junction2_chro_next = candidate[0][i+1].split(';')[1].split(':')[0]
                        #print(junction2_chro,junction1_chro_next)
                        if junction2_chro == 'gi|333031|lcl|HPV16REF.1|' and junction2_chro_next == 'gi|333031|lcl|HPV16REF.1|':

                            ax.plot([x-1,x-offset],[y,y],'-')
                            lines.append([(x-1, y),(x-offset,y)])
                        lc = mc.LineCollection(lines, linewidths=2,color = v_color)
                        ax.add_collection(lc)
                        lines = []
                        if junction2_chro != 'gi|333031|lcl|HPV16REF.1|' and junction2_chro_next != 'gi|333031|lcl|HPV16REF.1|':
                            ax.plot([x-1,x-offset],[y,y],'-')
                            lines.append([(x-1, y),(x-offset,y)])
                        lc = mc.LineCollection(lines, linewidths=2,color = h_color)
                        ax.add_collection(lc)
                        
    #                 print(i,candidate[3][i])
    #                 print(candidate[0])
                    #check looping node
                    copyNumberLoop = []
                    for eachLoop in candidate[3][i]:
                        k = 0

                        for eachNode in candidate[0][:i]:
                            if ';'.join(eachLoop.split(';')[0:2])  == eachNode:
                                #print(i,k,eachLoop,eachLoop.split(';')[1],candidate[0],candidate[1],candidate[3][i])
                                copyNumberLoop.append(str(i) + ';' + str(k))
                            k += 1
                    #plot the loops
                    countLoop = Counter(copyNumberLoop)
                    lines = []
                    for eachLoop in countLoop:
                        offset = 0.1
                        start = int(eachLoop.split(';')[0])
                        end = int(eachLoop.split(';')[1])
                        copyNumber = countLoop[eachLoop]
                        #color
                        if candidate[0][end].split(';')[1] == 'gi|333031|lcl|HPV16REF.1|' and candidate[0][start].split(';')[0] == 'gi|333031|lcl|HPV16REF.1|':
                            loopColor = v_color
                        else:
                            loopColor = h_color

                        lines.append([(lim-start-offset,y),(lim-start-offset,y+0.3)])
                        lines.append([(lim-start-offset,y+0.3),(lim-end-offset,y+0.3)])
                        lines.append([(lim-end-offset,y),(lim-end-offset,y+0.3)])
                        ax.text((lim-start-offset+lim-end-offset)/2, y+0.4, 'X'+str(copyNumber), zorder=100,
                            ha='center', va='center', weight='bold',
                            style='italic',size = 8)
                    if lines != []:
                        lc = mc.LineCollection(lines, linewidths=2,color = loopColor)
                        ax.add_collection(lc)


                    #print(Counter(copyNumberLoop))




                    i += 1
                j += 1


            ax.set_xlim(0,lim + 1)
            ax.set_ylim(0,lim + 1)
            fig.savefig(f'{outputPath}/{sampleName}.back.{m}.png')
            plt.close(fig)
        m += 1





