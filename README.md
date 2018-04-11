# A Dual Representation Simulated Annealing implementation

A C implementation of the Dual Representation Simulated Annealing (DRSA) [1] for the bandwidth minimization problem.
This implementation was made as close as possible to the algorithm description provided by [1].

## Compile and Run instructions 

* This code accept only graphs in the Matrix Market Format as input. See http://math.nist.gov/MatrixMarket for details.
* Go to the drsa folder and to compile type: $ make
* To run execute the exec_drsa file: $ ./exec_drsa

### Disclaimer

Note that I am not a DRSA author, so it is possible that this DRSA version has errors and/or discrepancies with the actual Torres-Jimenez et al. [1] DRSA algorithm. 

**[\[1\] Torres-Jimenez et al. A dual representation simulated annealing algorithm for the bandwidth minimization problem on graphs. Information Sciences 303 (2015) 33-49.](https://www.sciencedirect.com/science/article/pii/S0020025514011931)**
