#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:57:19 2021

@author: gs1511
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 12:55:45 2020

@author: gs1511
"""
class Node:
    def __init__(self, val):
        self.value = val
        self.leftChild = None
        self.rightChild = None
        self.indicator_val = None
        self.left_end = None
        self.right_end = None
        self.indicator = None # (1:refine, -1: coarsen, 0: do nothing)
        self.maxReached = 0 

    def insert(self, data):
        if self.value == data:
            return False
			
        elif self.value > data:
            if self.leftChild:
                return self.leftChild.insert(data)
            else:
                self.leftChild = Node(data)
                return True

        else:
            if self.rightChild:
                return self.rightChild.insert(data)
            else:
                self.rightChild = Node(data)
                return True 
    
    def find(self, data):
        if(self.value == data):
            return True
        elif self.value > data:
            if self.leftChild:
                return self.leftChild.find(data)
            else:
                return False
        else:
            if self.rightChild:
                return self.rightChild.find(data)
            else:
                return False
	
    def getSize(self):
        if self.leftChild and self.rightChild:
            return 1 + self.leftChild.getSize() + self.rightChild.getSize()
        elif self.leftChild:
            return 1 + self.leftChild.getSize()
        elif self.rightChild:
            return 1 + self.rightChild.getSize()
        else:
            return 1

    def getHeight(self):
        if self.leftChild and self.rightChild:
            return 1 + max(self.leftChild.getHeight(), self.rightChild.getHeight())
        elif self.leftChild:
            return 1 + self.leftChild.getHeight()
        elif self.rightChild:
            return 1 + self.rightChild.getHeight()
        else:
            return 1

    def preorder(self):
        if self:
            print (str(self.value))
            if self.leftChild:
                self.leftChild.preorder()
            if self.rightChild:
                self.rightChild.preorder()

    def postorder(self):
        if self:
            if self.leftChild:
                self.leftChild.postorder()
            if self.rightChild:
                self.rightChild.postorder()
            print (str(self.value))

    def inorder(self):
        if self:
            if self.leftChild:
                self.leftChild.inorder()
            print (str(self.value))
            if self.rightChild:
                self.rightChild.inorder()
    def refine(self):
        if self.indicator== 1 and self.maxRefReached == 0:
            if self.leftChild is None and self.rightChild is None:
                data_left = [self.value[0], .5*(self.value[0]+self.value[1])]
                data_right = [.5*(self.value[0]+self.value[1]), self.value[1]]
                self.leftChild = Node(data_left)
                self.rightChild = Node(data_right)
                return True
        return False
    
    def printLeaves(self):
        if self.leftChild is None and self.rightChild is None:
            print(str(self.value))
        else :
            self.leftChild.printLeaves()
            self.rightChild.printLeaves()
            
    
    def getLeaves(self,leaves_arr):
        #leaves_arr=[]
        if self.leftChild is None and self.rightChild is None:
            leaves_arr.append(self.value)
            #print(leaves_arr)
        else:
            self.leftChild.getLeaves(leaves_arr)
            self.rightChild.getLeaves(leaves_arr)
        #getLeaves(leaves_arr)
        #return leaves_arr
        
#def recursive_dfs(graph, node):
 #   seen = set()
# #   result = []

  #  def recursive_helper(node):
   #     for neighbor in graph[node]:
    #        if neighbor not in seen:
     #           result.append(neighbor)     # this line will be replaced below
      #          seen.add(neighbor)
       #         recursive_helper(neighbor)

    #recursive_helper(node)
    #return result
           
  
            
            

class Tree:
    def __init__(self):
        self.root = None

    def insert(self, data):
        if self.root:
            return self.root.insert(data)
        else:
            self.root = Node(data)
            return True
    

    def find(self, data):
        if self.root:
            return self.root.find(data)
        else:
            return False
	
    def getHeight(self):
        if self.root:
            return self.root.getHeight()
        else:
            return 0

    def getSize(self):
        if self.root:
            return self.root.getSize()
        else:
            return 0
	
    def remove(self, data):
        # empty tree
        if self.root is None:
            return False
        # data is in root node	
        elif self.root.value == data:
            if self.root.leftChild is None and self.root.rightChild is None:
                self.root = None
            elif self.root.leftChild and self.root.rightChild is None:
                self.root = self.root.leftChild
            elif self.root.leftChild is None and self.root.rightChild:
                self.root = self.root.rightChild
            elif self.root.leftChild and self.root.rightChild:
                delNodeParent = self.root
                delNode = self.root.rightChild
                while delNode.leftChild:
                    delNodeParent = delNode
                    delNode = delNode.leftChild
                
                if delNode.rightChild:
                    if delNodeParent.value > delNode.value:
                        delNodeParent.leftChild = delNode.rightChild
                    elif delNodeParent.value < delNode.value:
                        delNodeParent.rightChild = delNode.rightChild
                else:
                    if delNode.value < delNodeParent.value:
                        delNodeParent.leftChild = None
                    else:
                        delNodeParent.rightChild = None
                    self.root.value = delNode.value
						
            return True
		
        parent = None
        node = self.root
        
        # find node to remove
        while node and node.value != data:
            parent = node
            if data < node.value:
                node = node.leftChild
            elif data > node.value:
                node = node.rightChild
                
        # case 1: data not found
        if node is None or node.value != data:
            return False
		
        # case 2: remove-node has no children
        elif node.leftChild is None and node.rightChild is None:
            if data < parent.value:
                parent.leftChild = None
            else:
                parent.rightChild = None
            return True
		
        # case 3: remove-node has left child only
        elif node.leftChild and node.rightChild is None:
            if data < parent.value:
                parent.leftChild = node.leftChild
            else:
                parent.rightChild = node.leftChild
            return True
		
        # case 4: remove-node has right child only
        elif node.leftChild is None and node.rightChild:
            if data < parent.value:
                parent.leftChild = node.rightChild
            else:
                parent.rightChild = node.rightChild
            return True
		
        # case 5: remove-node has left and right children
        else:
            delNodeParent = node
            delNode = node.rightChild
            while delNode.leftChild:
                delNodeParent = delNode
                delNode = delNode.leftChild
			
            node.value = delNode.value
            if delNode.rightChild:
                if delNodeParent.value > delNode.value:
                    delNodeParent.leftChild = delNode.rightChild
                elif delNodeParent.value < delNode.value:
                   delNodeParent.rightChild = delNode.rightChild
            else:
                if delNode.value < delNodeParent.value:
                    delNodeParent.leftChild = None
                else:
                    delNodeParent.rightChild = None

    def preorder(self):
        if self.root is not None:
            print("PreOrder")
            self.root.preorder()
		
    def postorder(self):
        if self.root is not None:
            self.root.postorder()
	
    def inorder(self):
        if self.root is not None:
            print("InOrder")
            self.root.inorder()
    def printLeaves(self):
        if self.root is not None:
            self.root.printLeaves()
            
    def getLeaves(self):
        leaves = []
        if self.root is not None:
            self.root.getLeaves(leaves)
        return leaves
   
    
BST = Tree()
BST.insert(15)
BST.insert(8)
BST.insert(7)
BST.insert(9)
BST.insert(20)
BST.insert(17)
BST.insert(25)
BST.insert(21)
BST.insert(30)
BST.insert(1)
BST.insert(7.5)

BST.insert(.5)
BST.insert(1.5)
BST.insert(.25)
BST.insert(.75)
BST.insert(29)
BST.insert(31)






BST.inorder()
print('new')
BST.printLeaves()
m=[]
m=BST.getLeaves()
print(m)

BST2 = Tree()
BST2.insert(6)
BST2.insert(3)
BST2.insert(1)
BST2.insert(4)
BST2.insert(8)
BST2.insert(7)
BST2.insert(9)