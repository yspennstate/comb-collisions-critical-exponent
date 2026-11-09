# Infinitely many collisions on random comb graphs at the critical exponent

Yitzchak Shmalo

Paper: [paper/main.pdf](paper/main.pdf) (source [paper/main.tex](paper/main.tex)).

Let $\mathrm{Comb}(\mathbb Z,H)$ be the comb obtained by attaching to each $x\in\mathbb Z$ a tooth of height $H_x$,
with the $H_x$ i.i.d. and $P(H_x>z)\sim Cz^{-\gamma}$. De Ambroggio, Ng, Nitzschner and Scali (arXiv:2609.05343)
proved that two independent simple random walks on a typical such comb collide infinitely often almost surely when
$\gamma>1/3$ and only finitely often when $\gamma<1/3$, and conjectured infinitely many collisions at $\gamma=1/3$.
The paper proves the conjecture: at $\gamma=1/3$, for almost every comb, the two walks collide infinitely often
almost surely.

At the critical exponent the expected number of collisions in one scale is of order one and its second moment is of
order $\log N$. The Paley-Zygmund inequality then gives a collision in the $m$-th scale with conditional probability
at least $c/m$, and Levy's conditional Borel-Cantelli lemma applies because $\sum 1/m=\infty$.

## Code

No result of the paper depends on computation. `code/g_t_check.py` (Python 3, numpy and scipy) computes the
truncated backbone collision Green function $G_k(z)=\sum_{j\le k}\sum_y p_j(z,(y,0))^2$ of the lazy walk on
sampled combs with $\gamma=1/3$, for six starting points, as a numerical sanity check of the logarithmic growth;
`python g_t_check.py SEED LOG2T` appends to `code/g_t_check.out`, which holds the runs for seeds 1 and 4 up to
$k=2^{17}$.

## Funding

The research presented in this paper was supported by the European Research Council (ERC) under the European
Union's Horizon Europe research and innovation programme (grant agreement No. 101041711), by the Simons Foundation
as part of the Collaboration on the Mathematical and Scientific Foundations of Deep Learning, by Heights Labs, by
Convex Nexus Capital, by the Israel Science Foundation (grant number 2258/19), and by the Israel Science Foundation
(ISF Grant 4101/25). Computational resources were provided in part by Google Cloud through its education credits
program (credit 530550103).
