import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pandas import *
import pickle


Matrix_dict = pickle.load( open( "matrix.p", "rb" ) )

df = DataFrame(Matrix_dict).T.fillna(6)

dft = df.transpose()
f = dft.columns.values.tolist()
print f
dft.sort_values(by='33', inplace=True)

print dft

dft.to_csv('obsclusters.csv')

fig, ax = plt.subplots(1, 1)
cax = ax.imshow(dft, cmap='viridis', interpolation='nearest', origin='lower', aspect='auto')
ax.set_title('Heat map of observation clusters')
ax.set_xlabel("tile")
ax.set_ylabel("observations")
ax.legend()
plt.show()
