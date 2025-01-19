import numpy as np
from itertools import count

class sparse_l0:
    def __init__(self, A, b, penalty):
        self.A = A
        self.b = b
        self.penalty = penalty  # lambda
        self.L = np.max(np.abs(np.dot(A.T, A)))  # Lipschitz constant
        self.sample_size = A.shape[1]
        # initialization
        x = np.ones(self.sample_size)
        self.lengthx = len(x)
        self.x = x / np.sum(x)
        self.objective_list = [self.func_grad(variable=self.x)[0]]

    def ABPG_gain(self, tol=1e-6):
        rho = 1.2
        Gmin = 1e-2
        G = -1
        theta = 1
        gamma = 2
        z = self.x
        objective, a = self.func_grad(variable=self.x)
        for k in range(2000):
            Mk = max(G / rho, Gmin)
            theta_1 = theta
            for t in range(10000):
                G1 = Mk * rho ** t
                if k > 0:
                    theta = self.solve_theta(theta_1, gamma, G1 / G)
                y = (1 - theta) * self.x + theta * z
                fy, grad_fy = self.func_grad(variable=y)
                z1 = self.prox_map(z, grad_fy, G1 * theta ** (gamma - 1) * self.L)
                x1 = (1 - theta) * self.x + theta * z1
                fx, _ = self.func_grad(variable=x1)
                if fx <= fy + np.dot(grad_fy, (x1 - y)) + G1 * theta ** gamma * self.L * np.sum(z1 * np.log(z1 / z)):
                    break
            z = z1
            self.x = x1
            G = G1
            new_objective, a = self.func_grad(variable=x1)
            self.objective_list.append(new_objective)

            if k % 50 == 0:
                print(f"Iteration: {k}, objective value = {new_objective}, gradient norm = {np.linalg.norm(a)}.")
            if np.abs(new_objective - objective) < tol:
                print('The initialization step finishes!')
                print(f"The current iteration is {k}, the objective value is {new_objective}.")
                break
            objective = new_objective

    def step(self):
        print("Optimization begins! The initialization step starts. ")
        # Initialization step
        self.ABPG_gain()

        print('The L0BPG step starts.')
        # L0BPG step
        idx = 1
        objective = self.objective_value()
        for k in count():
            # BPG step
            y_new, _ = self.iteration()

            # sorting step
            idx, sample_size = self.search_dk(y_new, k, idx)
            remove_num = sample_size - idx
            if remove_num != 0:
                if k != 0:
                    # recover the vector with zero values
                    ynnew = np.zeros(self.lengthx)
                    # for i in range(len(ind)):
                    #     index = ind[i]
                    #     ynnew[index] = y_new[i]
                    ynnew[ind] = y_new
                    y_new = ynnew
                indx = np.argsort(y_new)
                ind = indx[self.lengthx - idx:]

                # removing step
                xm = y_new[ind]
                self.x = xm / np.sum(xm)
                self.A = self.A[:, ind]
            else:
                self.x = y_new
            new_objective = self.objective_value()
            self.objective_list.append(new_objective)
            if np.abs(objective - new_objective) < 1e-6:
                # recover final solution
                xnnew = np.zeros(self.lengthx)
                # for i in range(len(ind)):
                #     index = ind[i]
                #     xnnew[index] = self.x[i]
                xnnew[ind] = self.x
                self.x = xnnew
                print(f"Finish! The current iteration is {k}, the objective value is {new_objective} "
                      f"with nonzero number {np.sum(self.x>0)}.")
                break
            objective = new_objective
            if k % 10 == 0:
                print(f"Iteration: {k}, objective value = {new_objective}, nonzero number = {np.sum(self.x>0)}.")

    def func_grad(self, variable):
        '''
        calculate the objective value and the gradient
        :param variable:
        :return: objective value and gradient
        '''
        residual = np.dot(self.A, variable) - self.b
        grad = np.dot(self.A.T, residual)
        obj = 0.5 * np.linalg.norm(residual) ** 2
        return obj, grad

    def prox_map(self, variable, a, L):
        '''
        explicit expression for x
        '''
        x_without_n = variable * np.exp(-a / L)
        x = x_without_n / np.sum(x_without_n)
        return x

    def objective_value(self):
        obj = 0.5 * np.linalg.norm(np.dot(self.A, self.x) - self.b) ** 2 + self.penalty * np.sum(self.x != 0)
        return obj

    def solve_theta(self, theta, gamma, gainratio=1):
        # (1 - theta_k1) / theta_k1 ^ gamma = gainratio * 1 / theta_k ^ gamma
        ckg = theta ** gamma / gainratio
        cta = theta
        eps = 1e-6 * theta
        phi = cta ** gamma - ckg * (1 - cta)
        while abs(phi) > eps:
            drv = gamma * cta ** (gamma - 1) + ckg
            cta = cta - phi / drv
            phi = cta ** gamma - ckg * (1 - cta)
        return cta

    def iteration(self, num=1):
        y_new = self.x
        new_objective, a = self.func_grad(variable=self.x)
        for _ in range(num):
            y_new = self.prox_map(self.x, a, self.L)
            new_objective, a = self.func_grad(variable=y_new)
            # objective = new_objective
        return y_new, new_objective

    def search_dk(self, x_new, k, dk_last):
        '''
        search d_k in the sorting step
        '''
        a_sort = np.flip(np.sort(x_new[x_new > 0]))
        sample_size = len(a_sort)
        threshold = np.exp(self.penalty / self.L) - 1
        if k == 0:
            denominator = 0
            for m in range(1, sample_size):
                denominator += a_sort[m - 1]
                fraction = a_sort[m] / denominator
                if threshold > fraction:
                    return m, sample_size
        else:
            denominator = np.sum(a_sort)
            for m in range(dk_last, 0, -1):
                denominator -= a_sort[m - 1]
                fraction = a_sort[m - 1] / denominator
                if threshold <= fraction:
                    return m, sample_size
