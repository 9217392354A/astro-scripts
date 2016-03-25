#program for plotting functions
#chris Fuller, IRAM Thurs 27 June 2012

import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(1,20)
y = np.sqrt(x) / x


for i in range(1,13):
    ax = plt.subplot(6,2,i)    
    ax.plot(x,y*100.)
    ax.set_title(str(i))

plt.show() 