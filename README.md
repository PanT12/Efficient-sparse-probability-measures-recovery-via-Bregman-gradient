# Efficient sparse probability measures recovery via Bregman gradient
The complete paper can be accessed at <https://arxiv.org/pdf/2403.02861.pdf>. 

## Experiment Details
In the experiments, we consider synthetic data (Experiments 1-4) and real data (Experiments 5-6). 

### Experiment 1 (recovery accuracy)
In this experiment, we want to test the recovery accuracy of our algorithm. We maintain a constant signal-to-noise ratio (SNR) for data generation and compute the reconstruction SNR (RSNR) using the recovered solution.
The experiment demonstrates that our algorithm achieves high accuracy in a short time. Moreover, the solution proves to be robust, as evidenced by the low standard error. 

### Experiment 2 (support accuracy)
This experiment is conducted to determine if the solution derived from our algorithm can accurately recover the true support of the ground truth. To evaluate the performance of the algorithm, we compute various metrics, including accuracy, precision, recall, and the F1 score. Detailed definitions of these metrics are provided in the paper.
The experiment shows that our algorithm achieves high accuracy, precision, recall, and F1 value, indicating nearly perfect recovery. In addition, as the row dimension of matrix `A` increases, our algorithm demonstrates an increasingly remarkable capability to recover the ground truth vector. 

### Experiment 3 (efficacy of the L0BPG step)
The experiment indicates the efficacy of the L0BPG step in the algorithm. The algorithm would pick appropriate elements according to the sparsity penalty. Additionally, this experiment also highlights the minimal value control in the solution and the importance of achieving high accuracy during the initialization phase.

### Experiment 4 (Huber loss)
We consider the Huber loss as the loss function by introducing Salt and Pepper Impulse noise. We also compare it with the quadratic loss using our algorithm.

### Experiment 5 (hyperspectral unmixing)
In this experiment, we focus on a well-known region of the Cuprite dataset and try to recover the unknown abundance matrix. The experiment indicates that the abundance matrix recovered by our algorithm shows a high similarity to the Geological Reference Map. 

### Experiment 6 (portfolio optimization)
In this experiment, we consider portfolio optimization and use a benchmark dataset from the OR-Library. 

## Code Files Description
The data can be found under the file "data" and the result can be found under the file "result".

## Notes
1. The SUnSAL algorithm is provided by Jose Bioucas Dias at <http://www.lx.it.pt/~bioucas/code.htm>. 
2. Some codes in Experiment 5 are provided in <https://github.com/ricardoborsoi/MUA_SparseUnmixing/tree/master> [1].

## Reference
- [1] A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing. R.A. Borsoi, T. Imbiriba, J.C.M. Bermudez, C. Richard. IEEE Geoscience and Remote Sensing Letters, 2018.

