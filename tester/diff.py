import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

now = pd.read_csv('now.csv')
best = pd.read_csv('best.csv')

now_score = 0
best_score = 0

for i in range(now.shape[0]):
    if now['score'][i] < best['score'][i]:
        now_score += 1
    elif now['score'][i] > best['score'][i]:
        best_score += 1

print (now_score, best_score)
