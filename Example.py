import math
import numpy as np
from gurobipy import Model, GRB, LinExpr, quicksum
import time

m=Model("DW")
m.setParam('OutputFlag', 1) 
c=[14,8,11,7]
A=[[2.1,2.1,0.75,0.75],[0.5,0.5,0.5,0.5]]
b=[60,25]
A1=[[1,1,0,0],[1,0,0,0]]
A2=[[0,0,1,1],[0,0,1,0],[0,0,0,1]]
b1=[22,20]
b2=[12,15,25]
x=m.addVars(len(c),lb=0,vtype=GRB.CONTINUOUS,name="x")
m.setObjective(quicksum(c[i]*x[i] for i in range(len(c))),GRB.MAXIMIZE)
for i, row in enumerate(A):
    m.addConstr(quicksum(row[j] * x[j] for j in range(len(c))) <= b[i], name=f"c{i}")
for i, row in enumerate(A1):
    m.addConstr(quicksum(row[j] * x[j] for j in range(len(c))) <= b1[i], name=f"c1{i}")
for i, row in enumerate(A2):
    m.addConstr(quicksum(row[j] * x[j] for j in range(len(c))) <= b2[i], name=f"c2{i}")
t0=time.perf_counter()
m.optimize()
t1=time.perf_counter()
for v in m.getVars():
    print(f"{v.varName}: {v.x}")
print(f"Optimal Objective Value: {m.objVal}")
print(f"Total time: {t1-t0:.3f} seconds")
