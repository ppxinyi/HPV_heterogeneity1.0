import pysam
import re
import pandas as pd
from collections import Counter
import pickle
import os
from collections import defaultdict, deque
from typing import Dict, List, Tuple, Iterable, Optional
import pandas as pd
import re
#####extract junctions: extract junctions from nanopore bam file (minimap)
#bam: nanopore bam file
#chro: chromosome of virus integration
#ins1: integration start
#ins2: integration end
#virusChro: chromosome name of virus
#####
#return: junctionRes; format: key: readName, value:['readName','chro','junction','contigCoord','matchedLength','strand','direction'] ordered
#direction: l: left of matched part, r: right of matched part
# {'f1bf8c49-7d30-4d89-a257-d7a3bc37836a': [['f1bf8c49-7d30-4d89-a257-d7a3bc37836a',
#    '3',
#    189849673,
#    55190,
#    44661,
#    '-',
#    'l'],
def extract_junction(bam,chro,virusChro,ins1 = None,ins2 = None):
    sam = pysam.AlignmentFile(bam,'rb')
    count = 0
    SP = 0
    output = {}
    junctionRes = {}
    i = 0
    for read in sam.fetch(chro,ins1,ins2):
        i += 1
        if read.is_duplicate is True:
            continue
        if read.is_qcfail is True:
            continue
        if read.is_unmapped is True:
            continue
        if read.is_secondary is True:
            continue
        if read.mapq < 60:
            continue
        if read.has_tag('SA'):
            SP += 1
            if virusChro in read.get_tag('SA'):
            #if read.has_tag('SA'):
                junctionList = []
                #extract junction from primary cigar
                if read.is_reverse:
                    strand = '-'
                else:
                    strand = '+'
                if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                        print(read.infer_read_length())
                #calculate dulication
                matched_length = 0
                del_length = 0
                ins_length = 0
                for slot in read.cigar:
                    if slot[0] == 0: #matched length
                        matched_length += slot[1]
                    if slot[0] == 2: #deletion length
                        del_length += slot[1]
                    if slot[0] == 1:#insertion length
                        ins_length += slot[1]
                if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                    print(read.query_alignment_length,matched_length,del_length,ins_length)
                if read.cigar[0][0] == 0:
                    #add all matched, deleted basepairs 
                    humJunc = read.pos + matched_length + del_length
                    contigCoord = read.infer_read_length() - read.cigar[-1][1]
                    direction = 'r'
                    if strand == '-':
                        contigCoord = read.infer_read_length()-contigCoord
                        direction = 'l'
                    junctionList.append([read.qname,chro,humJunc,contigCoord,read.query_alignment_length,strand,direction])
                    if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                        print("from cigar", read.pos,humJunc)
                elif (read.cigar[0][0] == 4 or read.cigar[0][0] == 5) and (read.cigar[-1][0] == 4 or read.cigar[-1][0] == 5):
                    humJunc = read.pos
                    contigCoord = read.cigar[0][1]
                    direction1 = 'l'
                    humJunc2 = read.pos + matched_length + del_length
                    contigCoord2 = read.infer_read_length() - read.cigar[-1][1]
                    direction2 = 'r'
                    if strand == '-':
                        contigCoord = read.infer_read_length()-contigCoord
                        contigCoord2 = read.infer_read_length()-contigCoord2
                        direction1 = 'r'
                        direction2 = 'l'

                    junctionList.append([read.qname,chro,humJunc,contigCoord,read.query_alignment_length,strand,direction1])
                    junctionList.append([read.qname,chro,humJunc2,contigCoord2,read.query_alignment_length,strand,direction2])
                    if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                        print("from cigar",read.pos,humJunc,humJunc2)
                elif read.cigar[0][0] == 4 or read.cigar[0][0] == 5:
                    humJunc = read.pos
                    contigCoord = read.cigar[0][1]
                    direction = 'l'
                    if strand == '-':
                        contigCoord = read.infer_read_length()-contigCoord
                        direction = 'r'
                    junctionList.append([read.qname,chro,humJunc,contigCoord,read.query_alignment_length,strand,direction])
                    if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                        print("from cigar",read.pos,humJunc,humJunc2)
                #print(read.qname,read.cigar,read.get_tag('SA'))
                #extract junction from SA tag
                for slot in read.get_tag('SA').split(';')[:-1]:
                    
                    newChro = slot.split(',')[0]
                    start = int(slot.split(',')[1]) 
                    strand = slot.split(',')[2]
                    cigar = slot.split(',')[3]
                    nums = re.split('[A-Z]',cigar[:-1])
                    types = re.split('\d+',cigar)[1:]
                    #recode the matched part
                    j = 0
                    matched_length = 0
                    del_length = 0
                    nums = [int(x) for x in nums]
                    if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                        print(types,nums)
                    for each in types:
                        if each == 'M':
                            matched_length += nums[j]
                        if each == 'D':
                            del_length += nums[j]
                        j += 1
                    if types[0] == 'M': #start with matched part
                        rightjunction = start + matched_length + del_length
                        if rightjunction > 7906 and newChro == 'gi|333031|lcl|HPV16REF.1|':
                            rightjunction = rightjunction - 7906
                        rightCoord = read.infer_read_length() - nums[-1]
                        direction = 'r'
                        if strand == '-':
                            rightCoord = read.infer_read_length()-rightCoord
                            direction = 'l'
                        junctionList.append([read.qname,newChro,rightjunction,rightCoord,matched_length,strand,direction])
                        if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                            print('from SA',leftjunction,rightjunction)
                    elif (types[0] == 'H' or types[0] == 'S') and (types[-1] == 'H' or types[-1] == 'S'):
                        leftjunction = start
                        leftCoord = nums[0]
                        direction1 = 'l'
                        rightjunction = start + matched_length + del_length
                        if rightjunction > 7906 and newChro == 'gi|333031|lcl|HPV16REF.1|':
                            rightjunction = rightjunction - 7906
                        rightCoord = read.infer_read_length() - nums[-1]
                        direction2 = 'r'
                        if strand == '-':
                            rightCoord = read.infer_read_length()-rightCoord
                            leftCoord = read.infer_read_length()-leftCoord
                            direction1 = 'r'
                            direction2 = 'l'

                        junctionList.append([read.qname,newChro,leftjunction,leftCoord,matched_length,strand,direction1])
                        junctionList.append([read.qname,newChro,rightjunction,rightCoord,matched_length,strand,direction2])
                        if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                            print('from SA',leftjunction,rightjunction)
                    elif types[0] == 'H' or types[0] == 'S':
                        leftjunction = start
                        leftCoord = nums[0]
                        direction = 'l'
                        if strand == '-':
                            leftCoord = read.infer_read_length()-leftCoord
                            direction = 'r'
                        if read.qname == '436e9433-2402-4f77-9421-f203875eb89f':
                            print('from SA',leftjunction)
                        junctionList.append([read.qname,newChro,leftjunction,leftCoord,matched_length,strand,direction])
                    # for each_j in matched_j:
                    #     if each_j != 0:
                    #         leftjunction = start
                    #         leftCoord = sum(nums[0:each_j])
                    #         junctionList.append([read.qname,newChro,leftjunction,leftCoord,nums[each_j],strand,'l'])
                    #     if each_j != len(types)-1:
                    #         #add matched and deleted basepair
                    #         rightjunction = start + nums[each_j]
                    #         rightCoord = sum(nums[0:each_j+1])
                    #         if rightjunction > 7906 and newChro == 'gi|333031|lcl|HPV16REF.1|':
                    #             rightjunction = rightjunction - 7906
                    #         junctionList.append([read.qname,newChro,rightjunction,rightCoord,nums[each_j],strand,'r'])
                    #         if read.qname == '929f8266-ce0f-4c72-af79-0930909c3564':
                    #             print("from SA",rightjunction,slot,nums[each_j],each_j)
                if read.qname not in junctionRes:
                    junctionRes[read.qname] = junctionList                    
    #transfer to junction
    dataFrameList = []
    for read in junctionRes:
        dataFrame = pd.DataFrame(junctionRes[read],columns = ['readName','chro','junction','contigCoord','matchedLength','strand','direction'])
        dataFrameList.append(dataFrame)
        if read == '436e9433-2402-4f77-9421-f203875eb89f':
            print(junctionRes[read])

    return junctionRes,dataFrameList








#####merge junctions within numberBp, set numberBp
#junctionRes: output from extract_junction()
def merge_junction(junctionRes,numberBp=5):
    #transfer format of dictionary
    newJunctionRes = {}
    #print('junctionRes',junctionRes)
    for read in junctionRes:
        for each in junctionRes[read]:
            chro = each[1]
            if chro in newJunctionRes:
                newJunctionRes[chro].append(each)
            else:
                newJunctionRes[chro] = [each]
    
    #cluster integrations
    clusteredSiteRes = {}
    for eachChro in newJunctionRes:
        tempDf = pd.DataFrame(newJunctionRes[eachChro],columns = ['readName','chro','junction','contigCoord','matchedLength','strand','direction'])
        tempDf = tempDf.sort_values(by=['junction'],ignore_index = True)
#         k = 0
#         for each in tempDf['chro']:
#             if each == '3' and tempDf['junction'][k] == 189879024:
#                 print(tempDf.loc[[k]])
#             k += 1
        #Check circled HPV
        deltaList = [0]
#         if eachChro == 'gi|333031|lcl|HPV16REF.1|':
#             for i in range(0,len(tempDf)-1):
#                 #print(tempDf['junction'][i+1],tempDf['junction'][i],abs(tempDf['junction'][i+1]-tempDf['junction'][i])
#                 if abs(7906 - abs(tempDf['junction'][i+1]-tempDf['junction'][i])) > abs(tempDf['junction'][i+1]-tempDf['junction'][i]):
#                     deltaList.append(abs(tempDf['junction'][i+1]-tempDf['junction'][i]))
#                     print("circled:",tempDf['junction'][i+1],tempDf['junction'][i],abs(7906 - abs(tempDf['junction'][i+1]-tempDf['junction'][i])))
#                 else:
#                     deltaList.append(abs(7906 - abs(tempDf['junction'][i+1]-tempDf['junction'][i])))
#                     print("circled:",tempDf['junction'][i+1],tempDf['junction'][i],abs(7906 - abs(tempDf['junction'][i+1]-tempDf['junction'][i])))
                   
            
#         else:
        for i in range(0,len(tempDf)-1):
            deltaList.append(abs(tempDf['junction'][i+1]-tempDf['junction'][i]))
            
        cutList = []
        for i in range(len(deltaList)):
            if deltaList[i] > numberBp:
                cutList.append(i)
        #print(cutList)
        if cutList != []:
            clusteredSiteRes[eachChro] = [tempDf[0:cutList[0]].values.tolist(),tempDf[cutList[-1]:].values.tolist()]#add the first and last slots
            for i in range(len(cutList)-1):
                clusteredSiteRes[eachChro].append(tempDf[cutList[i]:cutList[i+1]].values.tolist())

        else:
            clusteredSiteRes[eachChro] = [tempDf.values.tolist()]
#     for chro in clusteredSiteRes:
#         i = 0
#         for eachCluster in clusteredSiteRes[chro]:
#             if chro == 'gi|333031|lcl|HPV16REF.1|':
#                 for read in eachCluster:
#                     if read[0] == '67195c32-12ea-4995-80fa-e0745c6cc641' and read[1] == 'gi|333031|lcl|HPV16REF.1|':
#                         print(read,i)
#             i += 1
            
    #replace the clustered junctionSite to be the same and reformat
   
    clusteredJunctionRes = {}
    for chro in clusteredSiteRes:
        for eachCluster in clusteredSiteRes[chro]:
            #replace clustered junctionSite to be the same
            if len(eachCluster) > 1:
                # if chro == '3':
                #     for slot in eachCluster:
                #         if slot[2] == 189889080:
                #             print(eachCluster)
                siteList = []
                #print(eachCluster)
                for each in eachCluster:
                    readName = each[0]
                    site = each[2]
                    siteList.append(site)
                countCluster = Counter(siteList)
                maxCount = max(countCluster.values())
                targetSite = 0
                for site in countCluster:
                    if countCluster[site] == maxCount:
                        targetSite = site
                        # if chro == '3':
                        #     print('targetSite',targetSite,chro)
                        break
                #write to the new dic
                
                for each in eachCluster:
                    readName = each[0]
                    if readName in clusteredJunctionRes:
                        # if chro == '3' and targetSite == 189889080:
                        #     print(each)
                        clusteredJunctionRes[readName].append(each[0:2] + [targetSite] + each[3:7] + [len(eachCluster)])
                    else:
                        clusteredJunctionRes[readName] = [each[0:2] + [targetSite] + each[3:7] + [len(eachCluster)]]
                    
            else:
                for each in eachCluster:
                    readName = each[0]
#                     if readName == 'fb410164-4e91-41b9-a2ad-c64df6a726cf':
#                         print(each)
                    if readName in clusteredJunctionRes:
                        clusteredJunctionRes[readName].append(each + [len(eachCluster)])
                    else:
                        clusteredJunctionRes[readName] = [each + [len(eachCluster)]]
    #print(clusteredJunctionRes)
    
#merge clustered junction HPV
    # mergedClusteredJunctionRes = {}

    # for read in clusteredJunctionRes:
    #     mergedClusteredJunctionRes[read] = []
    # #         if read == 'b1dbf8be-74a3-4e2c-b195-2465a8821676':
    #     #print(newJunctionRes[eachChro])
    #     tempDf = pd.DataFrame(clusteredJunctionRes[read],columns = ['readName','chro','junction','contigCoord','matchedLength','strand','direction','depth'])
    #     tempDf = tempDf.sort_values(by=['contigCoord'],ignore_index = True)
    #     i = 0
    #     #print(tempDf)
    #     for chro in tempDf['chro']:
    #         dup = False
    #         if i < len(tempDf)-1:
    #             if chro == 'gi|333031|lcl|HPV16REF.1|' and tempDf['chro'][i+1] == 'gi|333031|lcl|HPV16REF.1|':


    #                 #print(tempDf['junction'][i+1],tempDf['junction'][i],abs(tempDf['junction'][i+1]-tempDf['junction'][i])
    #                 if abs(7906 - abs(tempDf['junction'][i+1]-tempDf['junction'][i])) > abs(tempDf['junction'][i+1]-tempDf['junction'][i]):
    #                     delta = abs(tempDf['junction'][i+1]-tempDf['junction'][i])
    #                     #print("circled:",tempDf['junction'][i+1],tempDf['junction'][i],abs(7906 - abs(tempDf['junction'][i+1]-tempDf['junction'][i])))
    #                 else:
    #                     delta = abs(7906 - abs(tempDf['junction'][i+1]-tempDf['junction'][i]))
    #                 #print(delta)
    #                 if delta <= 5 and tempDf['contigCoord'][i+1]-tempDf['contigCoord'][i] <=5:
    #                     #print(tempDf)

    #                     dup = True
                        
    #             elif chro == tempDf['chro'][i+1]:
    #                 delta = abs(tempDf['junction'][i+1]-tempDf['junction'][i])
    #                 if delta <= 5 and tempDf['contigCoord'][i+1]-tempDf['contigCoord'][i] <=5:
    #                     dup = True
    #                     #print(tempDf)
    #                     #print("circled:",tempDf['junction'][i+1],tempDf['junction'][i],abs(7906 - abs(tempDf['junction'][i+1]-tempDf['junction'][i])))
    #             #merge junctions
    #         if dup:
    #             if tempDf['direction'][i] != tempDf['direction'][i+1]:
    #                 newLength = tempDf['matchedLength'][i] + tempDf['matchedLength'][i+1]
    #             else:
    #                 newLength = max(tempDf['matchedLength'][i],tempDf['matchedLength'][i+1])
    #             newJunction = [tempDf['readName'][i],tempDf['chro'][i],tempDf['junction'][i],tempDf['contigCoord'][i],\
    #                                newLength,\
    #                                tempDf['strand'][i],tempDf['direction'][i+1],max(tempDf['depth'][i],tempDf['depth'][i+1])]
                
    #             i += 1
    #         elif i < len(tempDf):
    #             newJunction = tempDf.iloc[i,:].tolist()
                
    #             mergedClusteredJunctionRes[read].append(newJunction)                       
                                   
    #         i += 1

    # for each in mergedClusteredJunctionRes:
    #     if each == 'b5517795-b0b6-4f68-8b92-7d07ec47b3be':
    #         # for junction in junctionListDic[each]:
    #         #     if '3:189887155' in junction:
    #         print(each,mergedClusteredJunctionRes[each])
    #         print('\n')
    # dataFrameList = []
    # for read in mergedClusteredJunctionRes:
    #     dataFrame = pd.DataFrame(mergedClusteredJunctionRes[read],columns = ['readName','chro','junction','contigCoord','matchedLength','strand','direction',
    #                                                                         'depth'])
    #     dataFrameList.append(dataFrame)
    dataFrameList = []
    for read in clusteredJunctionRes:
        dataFrame = pd.DataFrame(clusteredJunctionRes[read],columns = ['readName','chro','junction','contigCoord','matchedLength','strand','direction',
                                                                            'depth'])
        dataFrameList.append(dataFrame)
    junctionListDic = {}
    for read in dataFrameList:
        newRead = read.sort_values(by=['contigCoord','direction'],ascending = [True,False],ignore_index = True)
        if read['readName'][0] in junctionListDic:
            junctionListDic[read['readName'][0]].append(newRead)
        else:
            junctionListDic[read['readName'][0]] = [newRead]
    # for each in junctionListDic:
    #     for eachRead in junctionListDic[each]:
    #         if 3720 in list(eachRead['junction']):
    #             print(eachRead['junction'])
    #             print('\n')
    for each in junctionListDic:
        if each == '436e9433-2402-4f77-9421-f203875eb89f':
            # for junction in junctionListDic[each]:
            #     if '3:189887155' in junction:
            print('\n')
            print(each,junctionListDic[each])
            print('\n')
    #print(junctionListDic)
    
#     newClusteredJunctionRes = {}
#     for read in junctionListDic:
#         #if after merge next junction same, merge the next junctions
        
#         for each in junctionListDic[read]:
#             j = 0
#             while j < len(each):
                
# #                 if each['readName'][j] == '67195c32-12ea-4995-80fa-e0745c6cc641':
# #                     print(each.loc[[j]],j)
#                 if j + 1 < len(each.index):
#                     if each['chro'][j] == each['chro'][j+1] and each['junction'][j] == each['junction'][j+1]:
#                         m = 1 #skip all equal junction
#                         for k in range(j,len(each.index)-1):
#                             if each['chro'][k] == each['chro'][k+1] and each['junction'][k] == each['junction'][k+1]:
#                                 m += 1
#                                 print(k)
#                                 # newStart =each['contigCoord'][k]
#                                 # newLength = newLength + each['matchedLength'][k+1]
#                                 # newStrand = each['strand'][k+1]
#                                 # newDir = each['direction'][k+1]
#                                 # newDepth = each['depth'][k+1]
#                             else:
#                                 break
#                         j += m

#                 if j < len(each):
#                     newStart = each['contigCoord'][j]
#                     newLength = each['matchedLength'][j]
#                     newStrand = each['strand'][j]
#                     newDir = each['direction'][j]
#                     newDepth = each['depth'][j]
#                     if read in newClusteredJunctionRes:
#                         newClusteredJunctionRes[read].append(list(each.iloc[j,0:3]) + [newStart,newLength,newStrand,newDir,newDepth])
#                     else:
#                         newClusteredJunctionRes[read] = [list(each.iloc[j,0:3]) + [newStart,newLength,newStrand,newDir,newDepth]]
#                 j += 1

    
#     #change to dataframe 
#     dataFrameList = []
#     for read in newClusteredJunctionRes:
#         dataFrame = pd.DataFrame(newClusteredJunctionRes[read],columns = ['readName','chro','junction','contigCoord','matchedLength','strand','direction','depth'])
#         dataFrameList.append(dataFrame)
#     junctionListDic = {}
#     for read in dataFrameList:
#         newRead = read.sort_values(by=['contigCoord','direction'],ascending = [True,False],ignore_index = True)
#         if read['readName'][0] in junctionListDic:
#             junctionListDic[read['readName'][0]].append(newRead)
#         else:
#             junctionListDic[read['readName'][0]] = [newRead]
#     for each in junctionListDic:
#         if each == 'b5517795-b0b6-4f68-8b92-7d07ec47b3be':
#             # for junction in junctionListDic[each]:
#             #     if '3:189887155' in junction:
#             print(each,junctionListDic[each])
#             print('\n')
        
    return junctionListDic,clusteredJunctionRes


#####define junction
# #junctionListDic:output from merge_junction()
# def define_junction(junctionListDic):
#     #define junction
#     newJunctionListDic = {}
#     definedJunctionListDic = {}
#     for eachRead in junctionListDic:
#         i = 0
#         readList = junctionListDic[eachRead]
#         not_junction_index = []
#         for read in readList:
#             for eachReadName in read['readName']:
#             #filter continuous junction
#             #strategy: if matched length equal and direction is l + r, then continuous
#                 if i < len(read['readName'])-1:
#                     if read['direction'][i] == 'l' and  read['direction'][i+1] == 'r' and \
#                     read['matchedLength'][i+1] == read['matchedLength'][i]:
#                         junction_index = str(i) + ':' + str(i+1)
#                         not_junction_index.append(junction_index)
#                 i += 1
#             #not continuous = junction, need to devide into clean junction/ambigous junction
#             #clean junction: if contigCoord equal and direction is r + l, then clean
#             #ambigous junction: if contigCoord not equal but distance less than 100bp and direction is r + l, then ambigous
#             #undivided: else
#             i = 0
#             for each in read['readName']:
#                 if i < len(read['readName'])-1:
#                     junction_index = str(i) + ':' + str(i+1)
#                     newJunction = str(read['junction'][i]) + ':' + str(read['junction'][i+1])
#                     newChro = read['chro'][i] + ':' + read['chro'][i+1]
#                     newContigCoor = str(read['contigCoord'][i]) + ':' + str(read['contigCoord'][i+1])
#                     if (read['direction'][i]  !=  read['direction'][i+1]  )and \
#                     (read['contigCoord'][i+1] == read['contigCoord'][i]):
#                         newStrand = read['strand'][i] + ':' + read['strand'][i+1]
#                         if junction_index not in not_junction_index:
#                             if eachReadName in newJunctionListDic:
#                                 newJunctionListDic[eachReadName].append([eachReadName,newChro,newJunction,read['contigCoord'][i],newContigCoor,newStrand,'clean'])
#                             else:
#                                 newJunctionListDic[eachReadName] = [[eachReadName,newChro,newJunction,read['contigCoord'][i],newContigCoor,newStrand,'clean']]
#                     elif (read['direction'][i] !=  read['direction'][i+1])  and \
#                     ((read['contigCoord'][i+1] - read['contigCoord'][i]) < 100):
#                             newStrand = read['strand'][i] + ':' + read['strand'][i+1]
#                             #calibrate the order
#                             if i +2 < len(read['readName']):
#                                 if read['matchedLength'][i]==read['matchedLength'][i+2]:
#                                     newStrand = read['strand'][i+1] + ':' + read['strand'][i] 
#                                     newJunction = str(read['junction'][i+1]) + ':' + str(read['junction'][i])
#                                     newChro = read['chro'][i+1] + ':' + read['chro'][i]
#                             if i-1 > 0:
#                                 if read['matchedLength'][i+1] == read['matchedLength'][i-1]:
#                                     newStrand = read['strand'][i+1] + ':' + read['strand'][i] 
#                                     newJunction = str(read['junction'][i+1]) + ':' + str(read['junction'][i])
#                                     newChro = read['chro'][i+1] + ':' + read['chro'][i]
                           
                                
#                             if junction_index not in not_junction_index:
#                                 if eachReadName in newJunctionListDic:
#                                     newJunctionListDic[eachReadName].append([eachReadName,newChro,newJunction,read['contigCoord'][i],newContigCoor,newStrand,'ambiguous'])
#                                 else:
#                                     newJunctionListDic[eachReadName] = [[eachReadName,newChro,newJunction,read['contigCoord'][i],newContigCoor,newStrand,'ambiguous']]
#                     #else: undivided
#                 i += 1
        
#     #sort each read by contigCoord
#     dataFrameList = []
#     for read in newJunctionListDic:
#         dataFrame = pd.DataFrame(newJunctionListDic[read],columns = ['readName','chro','junction','contigCoord_1','contigCoord','strand','conf'])
#         dataFrameList.append(dataFrame)
#     definedJunctionListDic = {}
#     for read in dataFrameList:
#         newRead = read.sort_values(by=['contigCoord_1'],ascending = [True],ignore_index = True)
#         if read['readName'][0] in definedJunctionListDic:
#             definedJunctionListDic[read['readName'][0]].append(newRead)
#         else:
#             definedJunctionListDic[read['readName'][0]] = [newRead]
#     for each in definedJunctionListDic:
#         if each == '436e9433-2402-4f77-9421-f203875eb89f':
#             # for junction in junctionListDic[each]:
#             #     if '3:189887155' in junction:
#             print('\n')
#             print(each,definedJunctionListDic[each])
#             print('\n')

#     return definedJunctionListDic
###original function define avove 
###new define by Xinyi V2 add somenew
def define_junction(junctionListDic):
    """
    """
    newJunctionListDic = {}
    definedJunctionListDic = {}

    for read_key in junctionListDic:
        readList = junctionListDic[read_key]

        for read in readList:
            not_junction_index = []
            i = 0
            for _ in read['readName']:
                if i < len(read['readName']) - 1:
                    if (read['direction'][i] == 'l' and
                        read['direction'][i+1] == 'r' and
                        read['matchedLength'][i+1] == read['matchedLength'][i]):
                        not_junction_index.append(f"{i}:{i+1}")
                i += 1
            i = 0
            for _ in read['readName']:
                if i < len(read['readName']) - 1:
                    j = i + 1
                    junction_index = f"{i}:{j}"
                    chro1, chro2 = read['chro'][i], read['chro'][j]
                    junc1, junc2 = read['junction'][i], read['junction'][j]
                    contig1, contig2 = read['contigCoord'][i], read['contigCoord'][j]
                    strand1, strand2 = read['strand'][i], read['strand'][j]
                    dir1, dir2 = read['direction'][i], read['direction'][j]
                    ml1, ml2 = read['matchedLength'][i], read['matchedLength'][j]

                    newChro = f"{chro1}:{chro2}"
                    newJunction = f"{junc1}:{junc2}"
                    newContigCoor = f"{contig1}:{contig2}"
                    newStrand = f"{strand1}:{strand2}"
                    readName_out = read_key              
                    read_id_key_out = read_key           
                    matched_out = f"{ml1}:{ml2}"         

                    if (dir1 != dir2) and (contig2 == contig1):
                        if junction_index not in not_junction_index:
                            row = [
                                readName_out, newChro, newJunction,
                                contig1, newContigCoor, newStrand, 'clean',
                                read_id_key_out, matched_out
                            ]
                            newJunctionListDic.setdefault(readName_out, []).append(row)

                    elif (dir1 != dir2) and ((contig2 - contig1) < 100):
                        swapped = False
                        if (i + 2) < len(read['readName']):
                            if read['matchedLength'][i] == read['matchedLength'][i+2]:
                                (chro1, chro2)     = (chro2, chro1)
                                (junc1, junc2)     = (junc2, junc1)
                                (contig1, contig2) = (contig2, contig1)
                                (strand1, strand2) = (strand2, strand1)
                                (dir1, dir2)       = (dir2, dir1)
                                (ml1, ml2)         = (ml2, ml1)
                                swapped = True
                        if (i - 1) >= 0 and not swapped:
                            if read['matchedLength'][j] == read['matchedLength'][i-1]:
                                (chro1, chro2)     = (chro2, chro1)
                                (junc1, junc2)     = (junc2, junc1)
                                (contig1, contig2) = (contig2, contig1)
                                (strand1, strand2) = (strand2, strand1)
                                (dir1, dir2)       = (dir2, dir1)
                                (ml1, ml2)         = (ml2, ml1)

                        newChro       = f"{chro1}:{chro2}"
                        newJunction   = f"{junc1}:{junc2}"
                        newContigCoor = f"{contig1}:{contig2}"
                        newStrand     = f"{strand1}:{strand2}"
                        matched_out   = f"{ml1}:{ml2}"

                        if junction_index not in not_junction_index:
                            row = [
                                readName_out, newChro, newJunction,
                                contig1, newContigCoor, newStrand, 'ambiguous',
                                read_id_key_out, matched_out
                            ]
                            newJunctionListDic.setdefault(readName_out, []).append(row)
                i += 1

    cols = [
        'readName','chro','junction','contigCoord_1','contigCoord','strand','conf',
        'read_id_key','matchedLength'  # add new columns
    ]
    dataFrameList = []
    for readName_out in newJunctionListDic:
        df = pd.DataFrame(newJunctionListDic[readName_out], columns=cols)
        dataFrameList.append(df)

    definedJunctionListDic = {}
    for df in dataFrameList:
        newRead = df.sort_values(by=['contigCoord_1'], ascending=True, ignore_index=True)
        key = newRead['readName'][0]
        definedJunctionListDic.setdefault(key, []).append(newRead)

    return definedJunctionListDic
#33V2 only add one columns 
# def define_junction(junctionListDic):
#     newJunctionListDic = {}
#     definedJunctionListDic = {}

#     for eachRead in junctionListDic:
#         i = 0
#         readList = junctionListDic[eachRead]
#         not_junction_index = []

#         for read in readList:
#             for eachReadName in read['readName']:
#                 # filter continuous junction
#                 if i < len(read['readName']) - 1:
#                     if (read['direction'][i] == 'l' and
#                         read['direction'][i+1] == 'r' and
#                         read['matchedLength'][i+1] == read['matchedLength'][i]):
#                         junction_index = f"{i}:{i+1}"
#                         not_junction_index.append(junction_index)
#                 i += 1

#             # not continuous = junction
#             i = 0
#             for each in read['readName']:
#                 if i < len(read['readName']) - 1:
#                     junction_index = f"{i}:{i+1}"
#                     ml1, ml2 = read['matchedLength'][i], read['matchedLength'][i+1]
#                     matched_out = f"{ml1}:{ml2}"

#                     newJunction = f"{read['junction'][i]}:{read['junction'][i+1]}"
#                     newChro = f"{read['chro'][i]}:{read['chro'][i+1]}"
#                     newContigCoor = f"{read['contigCoord'][i]}:{read['contigCoord'][i+1]}"
#                     newStrand = f"{read['strand'][i]}:{read['strand'][i+1]}"

#                     # clean junction
#                     if (read['direction'][i] != read['direction'][i+1]) and \
#                        (read['contigCoord'][i+1] == read['contigCoord'][i]):
#                         if junction_index not in not_junction_index:
#                             row = [
#                                 eachReadName, newChro, newJunction,
#                                 read['contigCoord'][i], newContigCoor, newStrand,
#                                 'clean', matched_out
#                             ]
#                             newJunctionListDic.setdefault(eachReadName, []).append(row)

#                     # ambiguous junction
#                     elif (read['direction'][i] != read['direction'][i+1]) and \
#                          ((read['contigCoord'][i+1] - read['contigCoord'][i]) < 100):

#                         # calibrate the order
#                         if i + 2 < len(read['readName']):
#                             if read['matchedLength'][i] == read['matchedLength'][i+2]:
#                                 newStrand = f"{read['strand'][i+1]}:{read['strand'][i]}"
#                                 newJunction = f"{read['junction'][i+1]}:{read['junction'][i]}"
#                                 newChro = f"{read['chro'][i+1]}:{read['chro'][i]}"
#                                 matched_out = f"{ml2}:{ml1}"

#                         # only change here (>= 0)
#                         if (i - 1) >= 0:
#                             if read['matchedLength'][i+1] == read['matchedLength'][i-1]:
#                                 newStrand = f"{read['strand'][i+1]}:{read['strand'][i]}"
#                                 newJunction = f"{read['junction'][i+1]}:{read['junction'][i]}"
#                                 newChro = f"{read['chro'][i+1]}:{read['chro'][i]}"
#                                 matched_out = f"{ml2}:{ml1}"

#                         if junction_index not in not_junction_index:
#                             row = [
#                                 eachReadName, newChro, newJunction,
#                                 read['contigCoord'][i], newContigCoor, newStrand,
#                                 'ambiguous', matched_out
#                             ]
#                             newJunctionListDic.setdefault(eachReadName, []).append(row)
#                 i += 1

#     cols = [
#         'readName','chro','junction','contigCoord_1','contigCoord',
#         'strand','conf','matchedLength' 
#     ]
#     dataFrameList = []
#     for read in newJunctionListDic:
#         dataFrame = pd.DataFrame(newJunctionListDic[read], columns=cols)
#         dataFrameList.append(dataFrame)

#     definedJunctionListDic = {}
#     for read in dataFrameList:
#         newRead = read.sort_values(by=['contigCoord_1'], ascending=True, ignore_index=True)
#         key = newRead['readName'][0]
#         definedJunctionListDic.setdefault(key, []).append(newRead)

#     return definedJunctionListDic

###10.20Xinyi Add above 

###10.20look above

# #####count depth for defined junction
# #definedJunctionListDic: output from define_junction
# def cal_depth(definedJunctionListDic):
#     for read in definedJunctionListDic:
#         for df in definedJunctionListDic[read]:
#             df['junction']
        

#####output to intermedia file
#junctionListDic: clustered junction called from split reads
#outputPath: output directory
def write_junction(junctionListDic,sample,outputPath = './'):
    with open(outputPath + f'/{sample}.junctionListDic.pickle','wb') as output:
        pickle.dump(junctionListDic,output)

######read chromosome read coverage from samtools output
def read_cov(chro,file):
    with open(file) as inputFile:
        inputFile.readline()
        rows = inputFile.read().rsplit('\n')
        for row in rows:
            if row.split('\t')[0] == chro:
                cov = float(row.split('\t')[6])
    return cov


######add one function by XInyi

def flip_reads_by_nonhpv_minus(definedJunctionListDict):
    """
    """

    def is_hpv_token(token: str) -> bool:
        t = str(token).lower()
        return 'hpv' in t

    def row_has_nonhpv_minus(row) -> bool:
        chr_a, chr_b = map(str, row['chro'].split(':'))
        s_a, s_b = row['strand'].split(':')
        cond_a = (not is_hpv_token(chr_a)) and (s_a == '-')
        cond_b = (not is_hpv_token(chr_b)) and (s_b == '-')
        return cond_a or cond_b

    def swap_lr_fields(df: pd.DataFrame) -> pd.DataFrame:
        """(A:B → B:A)"""
        def swap_pair(x: str) -> str:
            if ':' in str(x):
                a, b = map(str, x.split(':'))
                return f"{b}:{a}"
            else:
                return x
        df = df.copy()
        for col in ['chro', 'junction', 'strand', 'contigCoord']:
            if col in df.columns:
                df[col] = df[col].apply(swap_pair)
        return df

    flipped_dict = {}

    for read_id, df_list in definedJunctionListDict.items():
        new_df_list = []
        for df in df_list:
            df = df.copy()
            df['flip_flag_row'] = df.apply(row_has_nonhpv_minus, axis=1)
            flip_flag = df['flip_flag_row'].any()

            df = df.sort_values('contigCoord_1', ascending=True, ignore_index=True)

            if flip_flag:
                df = df.iloc[::-1].reset_index(drop=True)
                df = swap_lr_fields(df)

            df = df.drop(columns=['flip_flag_row'])
            new_df_list.append(df)

        flipped_dict[read_id] = new_df_list

    return flipped_dict

####3V3 group read first and use longest read decide this set need reverse or not####

def flip_reads_by_weighted_host_strand_grouped(
    definedJunctionListDict: Dict[str, List[pd.DataFrame]],
    *,
    edge_overlap_threshold: int = 2,
    hpv_aliases: Iterable[str] = ('gi|333031|lcl|HPV16REF.1|',),
    drop_matchedLength_after_flip: bool = True,
) -> Dict[str, List[pd.DataFrame]]:
    """
    """
    def split_pair_to_str_ab(x: str) -> Tuple[str, str]:
        s = str(x)
        if ':' in s:
            a, b = s.split(':', 1)
            return a, b
        return s, s

    def split_pair_to_ints(x: str) -> Tuple[int, int]:
        a, b = split_pair_to_str_ab(x)
        def to_i(v):
            try:
                return int(float(str(v)))
            except Exception:
                return 0
        return to_i(a), to_i(b)

    POS_TOK = {'+', 'plus', 'pos', 'forward', 'fwd', '1', '+1'}
    NEG_TOK = {'-', 'minus', 'neg', 'reverse', 'rev', '-1'}
    def norm_strand_token(s: str) -> str:
        s = (str(s) or '').strip().lower()
        if s in POS_TOK: return '+'
        if s in NEG_TOK: return '-'
        if s in {'t', 'true', 'y', 'yes'}: return '+'
        if s in {'f', 'false', 'n', 'no'}: return '-'
        return '?'
    def norm_strand_pair(x: str) -> Tuple[str, str]:
        a, b = split_pair_to_str_ab(x)
        return norm_strand_token(a), norm_strand_token(b)
    
    alias_set = {x.lower() for x in (hpv_aliases or [])}
    def is_hpv_token(token: str) -> bool:
        s = str(token).strip().lower()
        return s in alias_set

    def row_has_hpv(row) -> bool:
        chr_a, chr_b = split_pair_to_str_ab(row.get('chro', ''))
        return is_hpv_token(chr_a) or is_hpv_token(chr_b)

    #reverse r and l :
    def swap_lr_fields(df: pd.DataFrame) -> pd.DataFrame:
        def swap_pair(x: str) -> str:
            a, b = split_pair_to_str_ab(x)
            return f"{b}:{a}"
        df = df.copy()
        for col in ['chro', 'junction', 'strand', 'contigCoord', 'matchedLength']:
            if col in df.columns:
                df[col] = df[col].astype(str).apply(swap_pair)
        return df

    def row_to_undirected_edge_key_stranded(row) -> str:
        chr_a, chr_b = split_pair_to_str_ab(row.get('chro', ''))
        jn_a, jn_b  = split_pair_to_str_ab(row.get('junction', ''))
        st_a, st_b  = norm_strand_pair(row.get('strand', '?:?'))
        va = f"{chr_a}@{jn_a}@{st_a}"
        vb = f"{chr_b}@{jn_b}@{st_b}"
        a_, b_ = sorted([va, vb])
        return f"{a_}|{b_}"
    ##create edge with HPV row
    read_edges = {}
    has_hpv_rows = {}
    for read_id, df_list in definedJunctionListDict.items():
        edges = set(); any_hpv = False
        for df in df_list or []:
            if df is None or len(df) == 0: continue
            if not {'chro', 'junction', 'strand'}.issubset(df.columns):
                continue
            mask = [row_has_hpv(row) for _, row in df.iterrows()]
            if not any(mask): continue
            any_hpv = True
            sub = df.loc[mask]
            for _, row in sub.iterrows():
                edges.add(row_to_undirected_edge_key_stranded(row))
        read_edges[read_id] = edges
        has_hpv_rows[read_id] = any_hpv

    hpv_reads = [r for r, v in has_hpv_rows.items() if v]
    adj = defaultdict(set)
    for r in hpv_reads:
        adj[r]  # ensure node

    for i in range(len(hpv_reads)):
        r1 = hpv_reads[i]; e1 = read_edges.get(r1, set())
        if not e1: continue
        for j in range(i + 1, len(hpv_reads)):
            r2 = hpv_reads[j]; e2 = read_edges.get(r2, set())
            if not e2: continue
            if len(e1 & e2) >= edge_overlap_threshold:
                adj[r1].add(r2); adj[r2].add(r1)

    components = []
    visit = set()
    for r in hpv_reads:
        if r in visit: continue
        comp = []
        q = deque([r]); visit.add(r)
        while q:
            cur = q.popleft(); comp.append(cur)
            for nb in adj[cur]:
                if nb not in visit:
                    visit.add(nb); q.append(nb)
        components.append(comp)

    def estimate_read_span(df_list):
        has_span = False; min_c1, max_c2 = None, None
        matched_sum = 0; edge_count = 0
        for df in df_list or []:
            if df is None or len(df) == 0: continue
            if 'chro' not in df.columns or 'junction' not in df.columns: continue
            mask = [row_has_hpv(row) for _, row in df.iterrows()]
            if not any(mask): continue
            sub = df.loc[mask]
            edge_count += len(sub)
            if 'contigCoord_1' in sub.columns and 'contigCoord_2' in sub.columns:
                c1 = pd.to_numeric(sub['contigCoord_1'], errors='coerce')
                c2 = pd.to_numeric(sub['contigCoord_2'], errors='coerce')
                if c1.notna().any() and c2.notna().any():
                    has_span = True
                    c1min, c2max = float(c1.min()), float(c2.max())
                    min_c1 = c1min if (min_c1 is None) else min(min_c1, c1min)
                    max_c2 = c2max if (max_c2 is None) else max(max_c2, c2max)
            if 'matchedLength' in sub.columns:
                for x in sub['matchedLength'].astype(str):
                    a, b = split_pair_to_ints(x)
                    matched_sum += (a + b)
        if has_span and min_c1 is not None and max_c2 is not None:
            return (max_c2 - min_c1 + 1.0), 'span'
        if matched_sum > 0:
            return float(matched_sum), 'matched_sum'
        return float(edge_count), 'edge_count'

    def compute_weighted_S(df_list):
        S_weight, has_weight, S_votes = 0, False, 0
        for df in df_list or []:
            if df is None or df.empty: continue
            if 'chro' not in df.columns or 'strand' not in df.columns: continue
            for _, row in df.iterrows():
                if not row_has_hpv(row): continue
                chr_a, chr_b = split_pair_to_str_ab(row.get('chro', ''))
                st_a, st_b = norm_strand_pair(row.get('strand', '?:?'))
                ml_a, ml_b = (0, 0)
                if 'matchedLength' in df.columns:
                    ml_a, ml_b = split_pair_to_ints(row.get('matchedLength', '0:0'))
                    if ml_a or ml_b: has_weight = True

                if not is_hpv_token(chr_a):
                    if has_weight:
                        if st_a == '+': S_weight += ml_a
                        elif st_a == '-': S_weight -= ml_a
                    else:
                        if st_a == '+': S_votes += 1
                        elif st_a == '-': S_votes -= 1

                if not is_hpv_token(chr_b):
                    if has_weight:
                        if st_b == '+': S_weight += ml_b
                        elif st_b == '-': S_weight -= ml_b
                    else:
                        if st_b == '+': S_votes += 1
                        elif st_b == '-': S_votes -= 1
        return S_weight if has_weight else S_votes

    flipped_dict = {
        r: [df.copy() if df is not None else df for df in (df_list or [])]
        for r, df_list in definedJunctionListDict.items()
    }

    flip_summary = [] 
    for comp_id, comp in enumerate(components):
        best_read, best_len, best_tag = None, float('-inf'), ''
        for r in comp:
            span, tag = estimate_read_span(definedJunctionListDict.get(r, []))
            if span > best_len:
                best_read, best_len, best_tag = r, span, tag
            elif span == best_len and best_read is not None:
                s_curr = abs(compute_weighted_S(definedJunctionListDict[r]))
                s_best = abs(compute_weighted_S(definedJunctionListDict[best_read]))
                if s_curr > s_best:
                    best_read, best_tag = r, tag

        need_flip = False
        if best_read is not None:
            S = compute_weighted_S(definedJunctionListDict[best_read])
            need_flip = (S < 0)

        if not need_flip:
            continue

        for r in comp:
            new_df_list = []
            for df in definedJunctionListDict.get(r, []) or []:
                if df is None or df.empty:
                    new_df_list.append(df); continue
                df2 = df.copy()
                if 'contigCoord_1' in df2.columns:
                    df2 = df2.sort_values('contigCoord_1', ascending=True, ignore_index=True)
                df2 = df2.iloc[::-1].reset_index(drop=True)
                df2 = swap_lr_fields(df2)
                if drop_matchedLength_after_flip and 'matchedLength' in df2.columns:
                    df2 = df2.drop(columns=['matchedLength'])
                new_df_list.append(df2)
            flipped_dict[r] = new_df_list

            flip_summary.append({
                "read_id": r,
                "component_id": comp_id,    
                "best_read": best_read,
                "S_value": S,
                "flip": bool(need_flip),
                "group_size": len(comp)
            })
    flip_summary_df = pd.DataFrame(flip_summary)
    return flipped_dict, flip_summary_df



#####check reverse Xinyi add


#### find path reverse
import pandas as pd
from typing import Dict, List
from collections import defaultdict

def detect_mirrors_after_flip(flipped_dict: Dict[str, List[pd.DataFrame]]):
    """
    """
    if not isinstance(flipped_dict, dict):
        raise TypeError("detect_mirrors_after_flip: expected a dict like {read_id: [df,...]}")

    def _flatten_junction_sequence(df_list: List[pd.DataFrame]) -> List[str]:
        seq = []
        for df in df_list or []:
            if df is None or df.empty or ('junction' not in df.columns):
                continue
            # if 'contigCoord_1' in df.columns:
            #     df = df.sort_values('contigCoord_1', ascending=True, ignore_index=True)
            seq.extend(df['junction'].astype(str).tolist())
        return seq

    def _canon_and_orientation(parts: List[str]):
        if not parts:
            return '', 'fwd'  
        fwd = tuple(parts)
        rev = tuple(reversed(parts))
        if fwd <= rev:
            return ' - '.join(fwd), 'fwd'
        else:
            return ' - '.join(rev), 'rev'

    rows = []
    for rid, dfs in flipped_dict.items():
        parts = _flatten_junction_sequence(dfs)
        canon_key, orient = _canon_and_orientation(parts)
        rows.append({
            'read_id': rid,
            'seq_len': len(parts),
            'canon_key': canon_key,
            'orientation_in_group': orient
        })

    df = pd.DataFrame(rows, columns=['read_id','seq_len','canon_key','orientation_in_group'])
    if df.empty or df['canon_key'].eq('').all():
        empty_summary = pd.DataFrame(columns=['canon_key','total_reads','n_fwd','n_rev'])
        empty_detail  = pd.DataFrame(columns=df.columns)
        return empty_summary, empty_detail

    stat = (df.groupby(['canon_key', 'orientation_in_group'])
              .size()
              .reset_index(name='n'))

    has_fwd = stat.loc[stat['orientation_in_group'].eq('fwd'), ['canon_key']].drop_duplicates()
    has_rev = stat.loc[stat['orientation_in_group'].eq('rev'), ['canon_key']].drop_duplicates()
    mirror_keys = has_fwd.merge(has_rev, on='canon_key', how='inner')

    mirrors_detail = (df.merge(mirror_keys, on='canon_key', how='inner')
                        .sort_values(['canon_key','orientation_in_group','read_id']))

    mirrors_summary = (mirrors_detail.groupby('canon_key', as_index=False)
                         .agg(total_reads=('read_id','count'),
                              n_fwd=('orientation_in_group', lambda s: (s=='fwd').sum()),
                              n_rev=('orientation_in_group', lambda s: (s=='rev').sum()))
                         .sort_values('total_reads', ascending=False))

    return mirrors_summary, mirrors_detail

