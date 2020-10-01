import pandas as pd

x = ["3", 4, 5]
y = sum([[value]*3 for value in x], [])

print(y)
