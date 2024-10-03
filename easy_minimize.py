import numpy as np
from scipy.optimize import minimize

def obj(r):
    sum = (r[0]-0)**2+(r[1]-1)**2+(r[2]-2)**2+(r[3]-3)**2
    return sum

eps = 1.0e-6

constrs = [
    {'type': 'eq', 'fun': lambda r: r[0]+r[1]+r[2]+r[3]-2}
]

res = minimize(obj, x0=np.ones(4), constraints=constrs)
print(res)