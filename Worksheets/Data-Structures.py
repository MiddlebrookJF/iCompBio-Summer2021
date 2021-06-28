import numpy as np
import pandas as pd

table = np.zeros((4, 2))
table[2,1] = 2

print(pd.DataFrame(table))
