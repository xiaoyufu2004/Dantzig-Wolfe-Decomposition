import gurobipy as gp
from gurobipy import GRB
import numpy as np
import time
# ==========================
# 数据定义
# ==========================
c1 = np.array([14, 8])     # 对应 x1, x2
c2 = np.array([11, 7])     # 对应 x3, x4

A1 = np.array([[2.1, 2.1],
               [0.5, 0.5]])
A2 = np.array([[0.75, 0.75],
               [0.5, 0.5]])

b = np.array([60, 25])

# 子问题局部约束
B1 = np.array([[1, 1],   # x1 + x2 <= 22
               [1, 0]])  # x1 <= 20
b1 = np.array([22, 20])

B2 = np.array([[1, 1],   # x3 + x4 <= 12
               [1, 0],   # x3 <= 15
               [0, 1]])  # x4 <= 25
b2 = np.array([12, 15, 25])

# ==========================
# 初始极点（至少给每块一个可行点）
# ==========================
X1_points = [np.array([0., 22.])]   # v1^1
X2_points = [np.array([0., 12.])]   # v2^1

# 初始极方向集合（开始为空）
X1_rays = []
X2_rays = []

# ==========================
# Master Problem（含极点和极方向）
# ==========================
def solve_master_problem(X1_pts, X2_pts, X1_rays, X2_rays, A1, A2, b, c1, c2):
    m = gp.Model("RMP")
    # λ for points; w for rays
    lam1 = m.addVars(len(X1_pts), lb=0, name="lam1")
    lam2 = m.addVars(len(X2_pts), lb=0, name="lam2")
    w1   = m.addVars(len(X1_rays), lb=0, name="w1")
    w2   = m.addVars(len(X2_rays), lb=0, name="w2")

    # 目标
    m.setObjective(
        gp.quicksum(c1 @ X1_pts[i]  * lam1[i] for i in range(len(X1_pts))) +
        gp.quicksum(c2 @ X2_pts[j]  * lam2[j] for j in range(len(X2_pts))) +
        gp.quicksum(c1 @ X1_rays[r] * w1[r]   for r in range(len(X1_rays))) +
        gp.quicksum(c2 @ X2_rays[s] * w2[s]   for s in range(len(X2_rays))),
        GRB.MAXIMIZE
    )

    # 耦合资源约束 A1*x1 + A2*x2 <= b
    for k in range(len(b)):
        m.addConstr(
            gp.quicksum(A1[k, :] @ X1_pts[i]  * lam1[i] for i in range(len(X1_pts))) +
            gp.quicksum(A2[k, :] @ X2_pts[j]  * lam2[j] for j in range(len(X2_pts))) +
            gp.quicksum(A1[k, :] @ X1_rays[r] * w1[r]   for r in range(len(X1_rays))) +
            gp.quicksum(A2[k, :] @ X2_rays[s] * w2[s]   for s in range(len(X2_rays)))
            <= b[k], name=f"res_{k}"
        )

    # 凸性约束：仅对 λ 生效
    m.addConstr(gp.quicksum(lam1[i] for i in range(len(X1_pts))) == 1, "conv1")
    m.addConstr(gp.quicksum(lam2[j] for j in range(len(X2_pts))) == 1, "conv2")

    m.setParam('OutputFlag', 0)
    m.setParam('Method', 1)  # 单纯形以便读取 .Pi
    m.optimize()

    if m.status != GRB.OPTIMAL:
        raise RuntimeError("RMP not optimal.")

    # 取对偶：耦合约束 π、凸性约束 μ
    pi  = [m.getConstrByName(f"res_{k}").Pi for k in range(len(b))]
    mu1 = m.getConstrByName("conv1").Pi
    mu2 = m.getConstrByName("conv2").Pi
    return (pi, mu1, mu2), m.objVal, m

# ==========================
# 极点定价（带 μ）
# ==========================
def solve_point_pricing(c, A, B, b, pi, mu):
    rc = c - A.T @ np.array(pi)
    m = gp.Model("PointPricing")
    x = m.addMVar(len(c), lb=0, name="x")
    # max (c - A^T pi)^T x - mu
    m.setObjective(rc @ x - mu, GRB.MAXIMIZE)
    m.addConstr(B @ x <= b)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status == GRB.OPTIMAL:
        return m.objVal, x.X
    return None, None

# ==========================
# 极方向定价（不含 μ；加归一化以取代表性方向）
# ==========================
def solve_ray_pricing(c, A, B, pi):
    rc = c - A.T @ np.array(pi)
    m = gp.Model("RayPricing")
    d = m.addMVar(len(c), lb=0, name="d")
    m.setObjective(rc @ d, GRB.MAXIMIZE)
    m.addConstr(B @ d <= 0)     # 齐次约束 => 方向
    m.addConstr(d.sum() <= 1)   # 归一化以避免无界，方便拿到一个 ray
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status == GRB.OPTIMAL:
        return m.objVal, d.X
    return None, None

# ==========================
# 迭代
# ==========================
tol = 1e-8
t0=time.perf_counter()
for it in range(50):
    (pi, mu1, mu2), obj, rmp = solve_master_problem(
        X1_points, X2_points, X1_rays, X2_rays, A1, A2, b, c1, c2
    )
    print(f"Iter {it+1:02d} | RMP obj={obj:.6f} | pi={pi} | mu=({mu1:.6f},{mu2:.6f})")

    # 定价：极点
    rc1_p, x1_new = solve_point_pricing(c1, A1, B1, b1, pi, mu1)
    rc2_p, x2_new = solve_point_pricing(c2, A2, B2, b2, pi, mu2)

    # 定价：极方向
    rc1_r, d1_new = solve_ray_pricing(c1, A1, B1, pi)
    rc2_r, d2_new = solve_ray_pricing(c2, A2, B2, pi)

    print(f"  RC point: X1={rc1_p:.6f}, X2={rc2_p:.6f} | RC ray: X1={rc1_r:.6f}, X2={rc2_r:.6f}")

    added = False
    # 只要 reduced cost > tol 就加入对应列
    if rc1_p is not None and rc1_p > tol:
        X1_points.append(x1_new); added = True
        print(f"  + add X1 point {x1_new}")
    if rc2_p is not None and rc2_p > tol:
        X2_points.append(x2_new); added = True
        print(f"  + add X2 point {x2_new}")
    if rc1_r is not None and rc1_r > tol:
        X1_rays.append(d1_new); added = True
        print(f"  + add X1 ray   {d1_new}")
    if rc2_r is not None and rc2_r > tol:
        X2_rays.append(d2_new); added = True
        print(f"  + add X2 ray   {d2_new}")

    if not added:
        print("No improving columns/rays. Optimal for DW RMP.\n")
        break
t1=time.perf_counter()
print(f"Total time: {t1-t0:.3f} seconds")
# （可选）从最终RMP还原 x1, x2：
# lam、w 的值可通过 rmp.getVarByName 读取并线性组合。
