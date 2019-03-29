import numpy as np
import matplotlib.pyplot as plt


# define column extraction 
def column(matrix, i):
        return [row[i] for row in matrix]


# read in file
solution = np.loadtxt('solution.plt', skiprows=1)


# split into variables
X = column(solution, 0)
u1 = column(solution, 1)
u2 = column(solution, 2)
u5 = column(solution, 3)
pres = column(solution, 4)


# plot stuff
plt.title('flow in a duct')
plt.plot(X,u5, '-')
plt.show()



