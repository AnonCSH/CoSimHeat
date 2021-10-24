# CoSimHeat：An Effective Heat Kernel Similarity Measure Based on Billion-Scale Network Topology  

## Tested Environment:

- Windows 10 
- Intel Core i7-6700 3.40GHz CPU and 64GB memory   
- visual studio 2019

## Datasets:

- DP：DBLP collaboration network
  - dblp.txt：graph
  - dictionary：the dictionary of id and name
  - DP_ground_truth.txt：ground truth
- EN：English word graph
  - EN_TS68_ground_truth.txt：ground truth
  - dictionary：the dictionary of id and word
  - edges.mtx：Normalized English graph (LFS)
  - outDegree.txt：the outDegree of graph
- EE：Email network of a EU institution
  - euAll.txt：graph
- WT：Wikipedia talk network
  - WT.txt：graph
- LJ：LiveJournal online social network
  - http://snap.stanford.edu/data/soc-LiveJournal1.html
- UK：A 2002 crawl of the .uk domain
  - https://sparse.tamu.edu/

Folders M and N under each dataset are randomly sampled nodes according to the output distribution

## Main:

(1) Accuracy:

```c++
// Search Quality on EN
acc_EN(c, theta, kmax); 
// Search Quality on DP
acc_DP(c, theta, kmax);
```

(2) CPU time :

```c++
timeTest(c, theta, kmax);
```

(3) Iterative Error:

```c++
errorTest(c, theta, kmax);
```
