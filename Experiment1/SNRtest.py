
from sparsel0 import *
from scipy.io import loadmat
from scipy.sparse import csr_matrix
import pandas as pd
import os
import timeit


if __name__ == '__main__':
    # A = loadmat('SNR_new/15/A.mat')['A']
    # b = loadmat('SNR_new/15/b1.mat')['b'].T[0]
    # x_true = loadmat('SNR_new/15/x1.mat')['x_true']
    # solution = sparse_l0(A, b, 1.5)
    # solution.step()


    numerator = 0
    denominator = 0
    numerator1 = 0
    denominator1 = 0
    SNR = []
    RSNR = []
    time = []
    num = 100

    for k in range(num):
        print(k)
        A = loadmat('Experiment1/50/A.mat')['A']
        b = loadmat(f'Experiment1/50/b{k + 1}.mat')['b'].T[0]
        x_true = loadmat(f'Experiment1/50/x{k + 1}.mat')['x_true']
        x_true = x_true.toarray().squeeze()

        SNR.append(10 * np.log10(np.linalg.norm(A @ x_true) ** 2 / np.linalg.norm(b - A @ x_true) ** 2))
        numerator += np.linalg.norm(A @ x_true) ** 2
        denominator += np.linalg.norm(b - A @ x_true) ** 2

        numerator1 += np.linalg.norm(x_true) ** 2

        lambda_ = 2

        # You need to implement or import the sparse_l0 function
        out = sparse_l0(A, b, lambda_)
        start_time = timeit.default_timer()
        out.step()
        end_time = timeit.default_timer()
        elapsed_time = end_time - start_time

        indt = np.where(x_true > 0)[0]
        indx = np.where(out.x > 0)[0]
        # common_elements = np.intersect1d(indt, indx)
        # print(common_elements)

        RSNR.append(10 * np.log10(np.linalg.norm(x_true) ** 2 / np.linalg.norm(out.x - x_true) ** 2))
        time.append(elapsed_time)
        denominator1 += np.linalg.norm(out.x - x_true) ** 2

    SNR_mean = 10 * np.log10(numerator / denominator)
    RSNR_mean = 10 * np.log10(numerator1 / denominator1)

    rowNames = [str(i + 1) for i in range(num)] + ['mean']
    result = np.column_stack((SNR, [SNR_mean] * num, RSNR, [RSNR_mean] * num, time, [2] * num))
    result = np.vstack((result, np.mean(result, axis=0)))
    colNames = ['SNR', 'SNR_mean', 'RSNR', 'RSNR_mean', 'time', 'lambda']
    result = pd.DataFrame(result, columns=colNames, index=rowNames)
    print(result)

