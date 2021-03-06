# Datasets

## Synthetic dataset 
The coefficients of polynomials are integers randomly chosen in the interval [-10,10], and the decision variable ![](https://latex.codecogs.com/gif.latex?V\in\{-1,1\}^n) with ![](https://latex.codecogs.com/svg.latex?2\leq%20n%20\leq%2010) for small-scale cases and ![](https://latex.codecogs.com/gif.latex?11\leq%20n\leq%2020) for relatively large-scale cases. The degree of polynomials d is chosen in the interval [2,6] which covers the most frequently used polynomials in real-world applications. 

## Benchmark dataset MQLib
This dataset https://github.com/niuyishuai/MQLib consists of 3296 large-scale, heterogeneous Max-Cut and Quadratic Unconstrained Binary Optimization (QUBO) problems, combining real-world problem instances and random problem instances from multiple random generators. Both of them are quadratic Boolean programs. More details can be found in Dunning, Gupta and Silberholz [1]

[1] I. Dunning, S. Gupta, and J. Silberholz, What works best when? a systematic evaluation of heuristics for max-cut and qubo, INFORMS Journal on Computing, 30 (2018), pp. 608--624
