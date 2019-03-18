import matplotlib.pyplot as plt
import numpy as np

def transform( data, N = None ):
    if ( N == None ):
        N = len(data)
        print(N)

    x = data[:N,0]
    y = data[:N,1]

    return np.log10(x), np.log10(y)

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

colors = ['blue', 'red', 'green', 'magenta', 'black']
orders = [4, 6, 8, 10, 12]

for order, color_ in zip(orders, colors):
    print('Order: {0}; color: {1}'.format(order, color_))
    filename = str(order) + 'order_diff_h.txt'

    data = np.loadtxt(filename)

    if ( order == 12 ):
        x, y = transform(data, 3)
        x2, y2 = transform(data, 20)
    else:
        x, y = transform( data, 20 )
    
    m, c = least_squares( x, y )
    print(m)
    xlin, ylin = approximate(m, c, x)

    plt.plot(xlin, ylin, linewidth = 0.5, color = color_) 
    plt.scatter(x, y, s = 5, color = color_)

    if ( order == 12 ):
        plt.scatter(x2, y2, s = 5, color = color_)

plt.show()

