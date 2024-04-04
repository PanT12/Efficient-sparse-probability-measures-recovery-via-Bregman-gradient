# Efficient sparse probability measures recovery via Bregman gradient
The complete paper can be accessed at <https://arxiv.org/pdf/2403.02861.pdf/>. 

## Experiment Details
In the experiments, we consider the least squares problem. 

### Experiment 1 (recovery accuracy)
In this experiment, we want to test the recovery accuracy of our algorithm. We maintain a constant signal-to-noise ratio (SNR) for data generation and compute the reconstruction SNR (RSNR) using the recovered solution.
The experiment demonstrates that our algorithm achieves high accuracy in a short time. Moreover, the solution proves to be robust, as evidenced by the low standard error. 

### Experiment 2 (support accuracy)
This experiment is conducted to determine if the solution derived from our algorithm can accurately recover the true support of the ground truth. To evaluate the performance of the algorithm, we compute various metrics, including accuracy, precision, recall, and the F1 score. Detailed definition of these metrics are provided in the paper.
The experiment shows that our algorithm achieves high accuracy, precision, recall and F1 value, indicating nearly perfect recovery. In addition, as the row dimension of matrix `A` increases, our algorithm demonstrates an increasingly remarkable capability to recover the ground truth vector. 

### Experiment 3 (efficacy of the L0BPG step)
The final experiment indicates the efficacy of the L0BPG step in the algorithm. The algorithm would pick appropriate elements according to the sparsity penalty. Additionally, this experiment also highlights the minimal value control in the solution and the importance of achieving high accuracy during the initialization phase.

## Code Files Description
The data and code can be found under the sections "Experiment 1", "Experiment 2", and "Experiment 3" files. The results of the algorithm can be viewed in the "Results" file.
