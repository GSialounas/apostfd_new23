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
            
    def get_leaf_nodes(self):
        leaves= []
        self._collect_leaf_nodes(self.root,leaves)
        return leaves
    
    def _collect_leaf_nodes(self,node,leaves):
        if node is not None:
            if (node.leftChild is None) and (node.rightChild is None):
                leaves.append(node.value)
            
            self._collect_leaf_nodes(node.leftChild,leaves)
            self._collect_leaf_nodes(node.rightChild,leaves)    

BST = Tree()
BST.insert(15)
print(BST.insert(8))
print(BST.insert(24))
print(BST.insert(5))
print(BST.insert(19))
print(BST.insert(30))
print(BST.insert(21))
print(BST.insert(9))
print(BST.insert(16))
print(BST.insert(32))
print(BST.insert(27))
BST.insert(26)
BST.remove(19)

BST.inorder()
