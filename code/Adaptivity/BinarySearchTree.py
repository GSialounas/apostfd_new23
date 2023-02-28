#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 21:48:31 2020

@author: gs1511
"""

class Node:
    def __init__(self,value=None):
        self.value = value
        self.left = None
        self.right = None
        self.children = []
    def add_child(self,object):
        self.children.append(object)
        
    
class BST:
    def __init__(self):
        self.root = None
    
    def insert(self, value):
        if self.root is None:
            self.root = Node(value)
        else:
            self._insert(value,self.root)
            
    def _insert(self,value, cur_node):
        print("cur_node value: " + str(cur_node.value) +" added value: "+str(value))

        if value < cur_node.value:
            print("got here")
            if cur_node.left is None:
                print("then here")
                cur_node.left = Node(value)
            else:
                print("or here")
                self._insert(value,cur_node.left)
        elif value>cur_node.value:
            if cur_node.right is None:
                cur_node.right = Node(value)
            else:
                self._insert(value, cur_node.right)
        else:
            print("No duplicates allowed in tree.")
    def remove(self,node):
        """We only remove leaf nodes and we do so for nodes of the same parent"""
        if node:
            if (node.left is None) and (node.right is None):
                node = None
            else:
                print("Node has children and cannot be removed")
                return False
            
    def find(self,value):
        if self.root:
            is_found = self._find(value,self.root)
            if is_found:
                return True
            return False
        else:
            return None
    def _find(self, value, cur_node):
        if value>cur_node.value and cur_node.right:
            return self._find(value,cur_node.right)
        elif value<cur_node.value and cur_node.left:
            return self._find(value,cur_node.left)
        if value ==cur_node.value:
            return True
    def get_leaf_nodes(self):
        leaves= []
        self._collect_leaf_nodes(self.root,leaves)
        return leaves
    
    def _collect_leaf_nodes(self,node,leaves):
        if node is not None:
            if (node.left is None) and (node.right is None):
                leaves.append(node.value)
            
            self._collect_leaf_nodes(node.left,leaves)
            self._collect_leaf_nodes(node.right,leaves)
    
        
    def print_tree(self, traversal_type):
        if traversal_type =="preorder":
            return self.preorder_print(self.root,"")
        elif traversal_type =="inorder":
            return self.inorder_print(self.root,"")   
        elif traversal_type =="postorder":
            return self.postorder_print(self.root,"")            
        else:
            print("U")
        
    def inorder_print(self, start, traversal):
        """Left->root->right"""
        if start:
            traversal = self.inorder_print(start.left,traversal)
            traversal +=(str(start.value) + "-")
            traversal = self.inorder_print(start.right,traversal)
        return traversal
    
    def postorder_print(self, start, traversal):
        """Left->right->root"""
        if start:
            traversal = self.inorder_print(start.left,traversal)
            traversal = self.inorder_print(start.right,traversal)
            traversal +=(str(start.value) + "-")
        return traversal
        
    def preorder_print(self, start,traversal):
        if start:
            traversal+=(str(start.value)+"-")
            traversal=self.preorder_print(start.left,traversal)
            traversal = self.preorder_print(start.right,traversal)
        return traversal
       

bst = BST()

bst.insert(5)
bst.insert(3)
bst.insert(7)
bst.insert(6)
bst.insert(8)     
bst.insert(4)
bst.insert(2)
#bst.insert(1)
#bst.insert(2.5)
#bst.insert(7.5)
#bst.insert(8.5)
#bst.insert(8.25)
#bst.insert(8.75)
#bst.insert(2.25)
#bst.insert(2.75)
#bst.insert(8.7)
#bst.insert(8.8)
print(bst.get_leaf_nodes())
bst.remove(bst.root.left.left)
print(bst.get_leaf_nodes())

     
                
            
        
            
        