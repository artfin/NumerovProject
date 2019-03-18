import matplotlib.pyplot as plt
import numpy as np

def transform( data, N ):
    x = data[1:,0][data[1:,0] < N]
    y = data[1:len(x)+1,1]

    return np.log(x), np.log(y)

def least_squares( x, y ):
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    return m, c

def approximate( m, c, x ):
    xmin = np.min(x)
    xmax = np.max(x)

    x = np.linspace(xmin, xmax, 100) 
    y = m*x + c
    return x, y

data = np.loadtxt("./8order_diff_k.txt")
x, y = transform( data, 10 )
m, c = least_squares( x, y )
print(m)
xlin, ylin = approximate(m, c, x)


#plt.text(1, 0, m, fontsize=12)
plt.plot(xlin, ylin, linewidth = 0.5, color = 'blue')
plt.scatter(x, y, s = 5, color = 'blue')
plt.show()
