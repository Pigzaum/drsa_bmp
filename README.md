# A Dual Representation Simulated Annealing implementation

A C implementation of the Dual Representation Simulated Annealing (DRSA) [[1](#references)] for the bandwidth minimization problem.
This implementation was made as close as possible to the algorithm description provided by [[1](#references)].

### Disclaimer

> Note that I am not a DRSA author, so it is possible that this DRSA version has errors and/or discrepancies with the actual Torres-Jimenez et al. [[1](#references)] DRSA algorithm. 

## Prerequisites

* GNU Make

* GCC

## Compile and Run instructions 

Go to the drsa folder and to compile type:

```sh
$ make
```

To run execute the lecm file:

```sh
$ ./exec_drsa -f <matrix file path>
```

This code accept only graphs in the Matrix Market Format as input. See http://math.nist.gov/MatrixMarket for details.


## References

**[\[1\] Torres-Jimenez et al. A dual representation simulated annealing algorithm for the bandwidth minimization problem on graphs. Information Sciences 303 (2015) 33-49.](https://www.sciencedirect.com/science/article/pii/S0020025514011931)**