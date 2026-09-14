# These classes do not need to be changed
    
# These classes do not need to be changed
    
class TreeNode:
    """ Class for holding nodes in the tree
    """
    def __init__(self, value=None):
        self._value = value
        self._children = []
        self._readNameList = {}
        self._endNode = []
        
    __all__ = ['value', 'children', 'add_child', 'add_readName','readNameList','endNode']
    
    #junction node
    @property
    def value(self):
        return self._value

    @property
    def children(self):
        return self._children
    
    #readName
    @property
    def readNameList(self):
        return self._readNameList
    
    #read duplication for end node
    @property
    def endNode(self):
        return self._endNode

    def add_child(self, child_node):
        self._children.append(child_node)
        
    def add_readName(self,readName):
        if self._value in self._readNameList:
            self._readNameList[self._value].append(readName)
        else:
            self._readNameList[self._value] = [readName]
            
    def add_endNode(self,endNode,readName):
        self._endNode.append(endNode + ';' + readName)
        
    
            
class BaseTree:
    """Main class for Tree
    This tree contains tree_node objects
    
    Private Attributes:
        nodes (list of tree_node): All of the nodes in the tree
    """
    def __init__(self):
        self.nodes = {} # Empty tree dictionary for references
        
    __all__ = ['add_node']
        
    def add_node(self, value, readName,parent_value=None,parent_readName=None):
        """ Adds a node to the tree
                by default there is no parent associated with a node
                
            Args:
                value (str): Value for the node (also this is the name of the node)
                parent_value (str): the name of the parent node
        """
        new_node = TreeNode(value)
        new_node.add_readName(readName)
        self.nodes[value] = new_node #Keep track of our objects in a dictionary
        
        
        if parent_value is not None:
            parent_node = self.nodes[parent_value]
            parent_node.add_child(new_node)
            
    def add_readName_in_node(self, value, readName,parent_value=None,parent_readName=None):
        """ Adds readname for existed node
                
                
            Args:
                value (str): Value for the node (also this is the name of the node)
                parent_value (str): the name of the parent node
        """
        if value == 'gi|333031|lcl|HPV16REF.1|:3720;3:189895061':
            print('\n',"readName",self.nodes)
        #search for this node in previous tree
        self.nodes[value].add_readName(readName)
        # if parent_value is not None:
        #     parent_node = self.nodes[parent_value]
        #     for child in parent_node.children:
        #         if child.value == value:
        #             child.add_readName(readName)
        #         # if value == 'gi|333031|lcl|HPV16REF.1|:1;gi|333031|lcl|HPV16REF.1|:1':
        #         #     print('\n',"readName",parent_value,child.readNameList,'\n')
                    
        # else:
        #     parent_node = self.nodes[value]
        #     parent_node.add_readName(readName)
        #     print('\n',"readName",parent_node.readNameList,'\n')
            
    def add_endNode_in_node(self, value, readName,parent_value=None,parent_readName=None):
        """ Adds end node for looping structure at the end node, e.g. ABCABCA
                
                
            Args:
                value (str): Value for the node (also this is the name of the node)
                parent_value (str): the name of the parent node
        """
        if parent_value is not None:
            parent_node = self.nodes[parent_value]
            parent_node.add_endNode(value,readName)
            
#     def add_depth(self,value,depth,parent_value=None,parent_readName=None):
#         """ Adds local read depth for existed node
                
                
#             Args:
#                 value (str): Value for the node (also this is the name of the node)
#                 parent_value (str): the name of the parent node
#         """
#         if parent_value is not None:
#             parent_node = self.nodes[parent_value]
#             for child in parent_node.children:
#                 if child.value == value:
#                     child.add_depth(depth)
#         else:
#             parent_node = self.nodes[value]
#             parent_node.add_depth(depth)
        
        
        
            

            
class Tree(BaseTree):
    def print_tree_depth_first(self, root_value):
        ''' This function will print our tree in a depth-first format.
        
        Args:
            root_value (str): The name of the root node in the tree
            depth: start of depth

        Returns:
            Prints tree depth-first       
        
        '''
        print(root_value)
        for child in self.nodes[root_value].children:
            
            self.print_tree_depth_first(child.value) 
            
            
    def print_tree_depth_first_read_name(self, root_value):
        ''' This function will print the read name list in our tree in a depth-first format.
        
        Args:
            root_value (str): The name of the root node in the tree

        Returns:
            Prints tree depth-first       
        
        '''
        print(root_value,self.nodes[root_value].readNameList[root_value])
        for child in self.nodes[root_value].children:
            #print(child.value,child.readNameList[child.value])
            self.print_tree_depth_first_read_name(child.value)
            
            
    def depth_first_search(self, root_value, search_value):
        ''' This function will perform a depth first search on our tree
        and will print every node value until it hits the item we are searching for.
        
        Args:
            root_value (str): The name of the root node in the tree
            search_value (str): The value that we are searching for

        Returns:
            found (bool): If the tree contains the value
            Prints tree depth-first up to search value
        
        
        '''
        print(root_value)
        if search_value == root_value:
            return True
            
        for child in self.nodes[root_value].children:
            if self.depth_first_search(child.value, search_value):
                return True
        
    def print_tree_breadth_first(self, root_value):
        ''' This function will print our tree in a breadth-first format.
        
        Args:
            root_value (str): The name of the root node to start printing from

        Returns:
            Prints tree breadth-first       
        
        '''
        queue = [root_value]
        queue2 = []
        i = 0
        while queue:
            current = queue.pop(0)
            print(current)
            #print(queue)
            if queue == []:
                i += 1
            print(i)
            for child in self.nodes[current].children:
                queue2.append(child.value)
            if queue == []:
                queue, queue2 = queue2, queue
                #print("child:",child.value,"parent:",current)
        #print(queue2)
    def breadth_first_search(self, root_value, search_value):
        ''' This function will perform a breadth first search on our tree
        and will print every node value until it hits the item we are searching for.
        
        Args:
            root_value (str): The name of the root node where the search will start
            search_value (str): The value that we are searching for

        Returns:
            found (bool): If the tree contains the value
            Prints tree depth-first up to search value
        
        
        '''
        queue = [root_value]
        
        while queue:
            current = queue.pop(0)
            print(current)
            if search_value == current:
                return True
            for child in self.nodes[current].children:
                queue.append(child.value)
        
        return False
    
#     def binaryTreePaths(self,root_value):
        #resList = []
    def walk_assembly(self,root_value,res,score,readNames,endNodeList,resList):
        ''' This function will perform a walk to rank candidates of possible structures

        Args:
            root_value (str): The name of the root node where the search will start

        Returns:
            Lists of nodes
            Prints candidate of assemblies, including junctions and reads


        '''


        #store number of reads and number of segements in a dataframe
        res.append(root_value)
        reads = self.nodes[root_value].readNameList[root_value]
        readNameList = []
        for each in reads:
            readNameList.append(each.split(';')[0])
        numReads = len(set(readNameList))
        numSegs = len(readNameList)
        score.append([numReads,numSegs])
        readNames.append(self.nodes[root_value].readNameList[root_value])
        endNode = self.nodes[root_value].endNode
        endNodeList.append(endNode)
        
        if self.nodes[root_value].children == []:
#             print(res)
#             print(score)
            resList.append([res.copy(),score.copy(),readNames.copy(),endNodeList.copy()])
            #print(resList)
        for child in self.nodes[root_value].children:
            #res.append(child.value)
            self.walk_assembly(child.value,res,score,readNames,endNodeList,resList)
            res.pop()
            score.pop()
            endNodeList.pop()
            readNames.pop()
            