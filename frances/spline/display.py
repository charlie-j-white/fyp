import numpy as np
import matplotlib.pyplot as plt


# define column extraction 
def column(matrix, i):
        return [row[i] for row in matrix]


# read in file
solution = np.loadtxt('transect.plt', skiprows=1)


# split into variables
X = column(solution, 0)
u1 = column(solution, 1)


# plot stuff
plt.title('spline')
a, = plt.plot(X,u1, '.',label='x,y')
plt.legend(handles=[a])
plt.show()



