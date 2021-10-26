# CoSimHeat：An Effective Heat Kernel Similarity Measure Based on Billion-Scale Network Topology  
## Environment:

- Windows 10 
- Intel Core i7-6700 3.40GHz CPU and 64GB memory   
- Visual Studio C++ 2019

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

## Compared Algorithms:

We compared a series of algorithms, whose implementation files are shown below:

- CoSimHeat (our) : CoSimHeat.h 
- CoSimRank : CoSR.h
- RP-CoSim : RPCoSim.h
- SimRank : ExactSim.h
- SimRank* : SRS.h
- ASCOS : ASCOS.h

## Datasets:

- DP：DBLP collaboration network
  - dblp.txt：graph
  - dictionary：the dictionary  (id, name)
  - DP_ground_truth.txt：ground truth
- EN：English word graph
  - EN_TS68_ground_truth.txt：ground truth
  - dictionary：the dictionary  (id, word)
  - edges.mtx：Normalized English graph
  - outDegree.txt：the outDegree of each node on graph
- EE：Email network of a EU institution
  - euAll.txt：graph
- WT：Wikipedia talk network
  - WT.txt：graph
- LJ：LiveJournal online social network
  - http://snap.stanford.edu/data/soc-LiveJournal1.html
- UK：A 2002 crawl of the .uk domain
  - https://suitesparse-collection-website.herokuapp.com/MM/LAW/uk-2002.tar.gz

Folders M and N under each dataset store the query nodes that are randomly sampled according to the degree distribution of each graph.
