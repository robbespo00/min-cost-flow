# **Computational Mathematics for Learning and Data Analysis**

<p align="center">
  <img height=500px width=550px src="https://user-images.githubusercontent.com/73891662/148542540-b0cb8277-d883-46a2-8de1-3058101c0372.PNG">
</p>

# **Getting Started**
The following repository contains the non ML project 27 of the CM course.<br />

## Prerequisites
Given N the number of nodes and E the number of edges of the graph, the instances we use come from a `.txt` file and have the following structure: the first N rows are formatted like that: `'n' | node | supply`, while the rows from the (N+1)-th to the last (E+N) have this format: `'a' | source_node | destinatation_node | cost`.<br />

## Execution
The main function is `main(filename, mode, generate, distribution)` where:
* `filename` is the path of the file containing the instance;
* `mode` is the execution mode ('residual', 'rate', 'minres', 'precond');
* `generate` is a boolean variable which indicates if the function has to generate or to load the vector D;
* `distribution` is the probability distribution of data in the vector D choosen between 7 fixed distributions.
