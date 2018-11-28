# Minimum_Spanning_Tree

Author          : Yingqiao Zheng
Created         : September 28, 2018
Last Modified   : October 2, 2018

Affiliation     : Georgia Institute of Technology


Description
-------------

This program is to compute Minimum Spanning Tree (MST) for a specific graph using DFS, and by reading an extra file of new edges, recompute MST for the graph (without perform DFS to the whole again), and finally output the total MST weight to file after processing each new edge.


Installation
------------

To install, simply run

    g++ RunExperiment.c -o MST -std=c++11


Execution
-----------

Assuming your executable is called "MST", run it using

    ./MST <original_graph_path> <extra_edge_file_path> <output_file_path>


