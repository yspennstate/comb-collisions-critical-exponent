"""Truncated backbone collision Green function on a random comb with gamma = 1/3.

For the lazy simple random walk on Comb(Z, H), with i.i.d. teeth P(H >= h) = h^{-gamma}, this computes
G_k(z) = sum_{j<=k} sum_y p_j(z,(y,0))^2 at k = 2^4, ..., T for six starting points z. At the critical
exponent G_k is expected to grow like log k, i.e. by a roughly constant amount per doubling of k.
Sanity check only; no result of the paper depends on it.

Usage: python g_t_check.py SEED LOG2T   (appends to g_t_check.out next to this file)"""
import os
import sys
import numpy as np
import scipy.sparse as sp

seed, l2t = int(sys.argv[1]), int(sys.argv[2])
g = 1.0 / 3.0
T = 2 ** l2t
rng = np.random.default_rng(seed)
rangeT = T ** ((g + 1) / 4)
L = int(6 * rangeT) + 10
Hmax = int(4 * T ** 0.5) + 10
zs = np.arange(-L, L + 1)
H = np.minimum(np.floor(rng.random(len(zs)) ** (-1.0 / g)).astype(np.int64), Hmax)
off = np.zeros(len(zs) + 1, dtype=np.int64)
off[1:] = np.cumsum(H + 1)
n = off[-1]
rows, cols, vals = [], [], []


def add(i, nbrs):
    w = 0.5 / len(nbrs)
    rows.append(i); cols.append(i); vals.append(0.5)
    for j in nbrs:
        rows.append(i); cols.append(j); vals.append(w)


for a in range(len(zs)):
    b = off[a]
    nb = []
    if a > 0: nb.append(off[a - 1])
    if a < len(zs) - 1: nb.append(off[a + 1])
    if H[a] >= 1: nb.append(b + 1)
    add(b, nb)
    for y in range(1, H[a] + 1):
        add(b + y, [b + y - 1] + ([b + y + 1] if y < H[a] else []))
P = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
PT = P.T.tocsr()
back = off[:-1]
# starts: origin, the base of the tallest tooth within rangeT, and 4 random backbone sites within rangeT
near = np.where(np.abs(zs) <= rangeT)[0]
starts = [L, near[np.argmax(H[near])]] + list(rng.choice(near, 4, replace=False))
X = np.zeros((n, len(starts)))
for c, a in enumerate(starts):
    X[off[a], c] = 1.0
G = np.zeros(len(starts))
out = []
for k in range(1, T + 1):
    X = PT @ X
    G += (X[back] ** 2).sum(0)
    if k & (k - 1) == 0 and k >= 16:
        out.append("seed %d k 2^%d G %s" % (seed, k.bit_length() - 1, " ".join("%.4f" % v for v in G)))
hdr = "seed %d T 2^%d states %d L %d Hmax %d start_z %s start_H %s" % (
    seed, l2t, n, L, Hmax, [int(zs[a]) for a in starts], [int(H[a]) for a in starts])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "g_t_check.out"), "a") as f:
    f.write(hdr + "\n" + "\n".join(out) + "\n")
print(hdr)
