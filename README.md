# GAP code for "Residual Subnormal Defect Two: Affine Criteria and Extremal Fitting Height"

This repository contains the [GAP](https://www.gap-system.org) programs used for the computations
reported in

> L. Gong, Y. Yan and E. Cosme Llópez, *Residual Subnormal Defect Two: Affine Criteria and Extremal Fitting Height*.

together with the output of the runs (`logs/`).

For a finite group $G$, the tested residuals are the nilpotent residuals $H^{\mathcal N}$ of the subnormal
subgroups $H$ of $G$. The residual defect $b_{\mathcal N}(G)$ is the largest subnormal defect of a tested
residual, $b_{\mathrm{sn}}(G)$ is the largest subnormal defect of a subnormal subgroup, $\mathrm W(G)$ is the
intersection of the normalizers of the tested residuals, $\ell_{\mathrm W}(G)$ is the length of the upper
$\mathrm W$-series and $h(G)$ is the Fitting height. Notation follows the paper, and theorem and section
numbers refer to it. References of the form [Part I, ...] are to the companion paper
L. Gong, Y. Yan and E. Cosme Llópez, *Wielandt subgroups of nilpotent residuals and Fitting-height bounds*,
[arXiv:2610.28082](https://arxiv.org/abs/2610.28082).

The paper itself relies on a computer search only for Example 3.11. The other scripts check
statements of the paper on small cases, as an independent test of the proofs.

## Requirements

GAP 4.16 or later, with the SmallGroups library (included in the standard GAP distribution).
The runs in `logs/` were made with GAP 4.16. The shell script needs `bash` and `xargs`. All commands are run
from the root of the repository.

## Files

| File | Contents |
|---|---|
| `gap/lib.g` | Common functions: nilpotent residual, subnormal defect, subnormal subgroups up to conjugacy (`SubnormalReps`), $O^p(A)$ (`PResidual`), $b_{\mathcal N}(G)$ and $\mathrm W(G)$ through $G/\mathrm F(G)$ (`ResidualData`) or by brute force, also giving $b_{\mathrm{sn}}(G)$ (`ResidualDataBrute`), $\ell_{\mathrm W}(G)$, Fitting height, and affine groups $V\rtimes M$ as pc groups (`AffinePcGroup`). Must be loaded before the other scripts. |
| `gap/smallest.g` | Example 3.11: every group of order less than 72 has $b_{\mathcal N}\le1$; $S_3\wr C_2$ = `SmallGroup(72,40)` is the only group of order 72 with $b_{\mathcal N}=2$, and it has $b_{\mathrm{sn}}=2$; in the affine model, $\mathrm W(G)=V\rtimes Z(M)$. |
| `gap/affine.g` | Section 3: Theorem 3.3 (the residual formula $H^{\mathcal N}=[V,O^p(A)]A^{\mathcal N}$), Theorem 3.1 (the criterion for $b_{\mathcal N}\le2$), Theorem 3.5 (the formula for $\mathrm W(G)$) and Proposition 3.8, for eleven affine groups. When $\dim V\le4$ the formula is checked for all subnormal subgroups by brute force; for larger $V$, where every subspace of $V$ is subnormal, it is checked for the subgroups $VA$. |
| `gap/uniform.g` | Theorem 3.10: the condition on $M$ alone against $b_{\mathcal N}(\mathbb F_p[M]\rtimes M)$, for the regular modules of seven pairs $(M,p)$. |
| `gap/blocks.g` | Lemmas 4.3, 4.4, 5.4 and 5.5: the blocks $\mathrm{AGL}(2,3)$, $B_{13}$, $B_{17}$ and $\mathbb F_{29}^3\rtimes(C_7\rtimes C_3)$. |
| `gap/construction.g` | Section 5: three small instances of the inverse-image construction of Theorem 5.10, checking its conclusions and those of Proposition 5.8. |
| `gap/obstruction.g` | Section 6: Theorem 6.1 on $E=S_4\wr C_2$, $p=3$, including a tested residual of defect three, and the hypotheses of Corollary 6.2. |
| `gap/bryce.g` | Theorem 7.10: the explicit group $X_2$ of Bryce's tower, of order $75\cdot17^{12}$, with $b_{\mathcal N}=2$ and a subnormal subgroup of defect three. |
| `run_checks.sh` | Reproduces all the logs. |
| `logs/` | Output of the runs. |

## Usage

```
gap -q gap/lib.g gap/smallest.g
gap -q -c 'long:=true;;' gap/lib.g gap/construction.g   # includes the long case B
```

To reproduce all the logs (4 parallel processes, about an hour and a half, mostly Case B of `construction.g`):

```
./run_checks.sh
```

Each script prints one line per checked statement or group, with the computed values and `true` or
`false` for each claim, followed by `FINAL`. Where a value is expected, the line says so.

## Results

All the claims checked hold.

- **Example 3.11** (`smallest.g`): among the groups of order at most 72, only `SmallGroup(72,40)` has
  $b_{\mathcal N}\ge2$; it has $b_{\mathcal N}=b_{\mathrm{sn}}=2$, and $|G:\mathrm W(G)|=4$.
- **Section 3** (`affine.g`, `uniform.g`): Theorems 3.1, 3.3, 3.5 and 3.10 and Proposition 3.8 hold in
  all the cases, which include groups with $b_{\mathcal N}=1$ and $b_{\mathcal N}=2$, and groups with
  $b_{\mathcal N}(V\rtimes M)=b_{\mathcal N}(M)+1$. In all the pairs of `uniform.g` the condition of
  Theorem 3.10 holds; a case where it fails is in `obstruction.g`.
- **Blocks** (`blocks.g`): both blocks are $\mathrm W$-groups of the stated height, $K^{\mathcal N}=P$, and
  no proper subnormal subgroup $L$ of $K$ satisfies $LP=K$. For $\ell=13,17$ the single-coordinate
  vectors of $V_\ell^*$ form $(\ell^2-1)/8$ regular orbits of cyclic vectors, and for the odd block
  they form $(t-1)/r=4$ such orbits.
- **Section 5** (`construction.g`): in the three cases, $|G|=|K|^{|E|}|D||E|$, $h(G)=h(E)+h(B)$,
  $\mathrm F(G)=V^{\Omega}$, $b_{\mathcal N}(G)=2$, $N\le\mathrm W(G)$ and
  $\ell_{\mathrm W}(G)=\ell_{\mathrm W}(E)+1$.
- **Section 6** (`obstruction.g`): for $E=S_4\wr C_2$ and $p=3$ the condition of Theorem 3.10 fails, and
  already the permutation module $\mathbb F_3^8$ gives a tested residual of defect three. The
  hypotheses of Corollary 6.2 hold for the three blocks.
- **Theorem 7.10** (`bryce.g`): $|G|=75\cdot17^{12}$, $\mathrm F(G)=V$, $h(G)=3$, $\ell_{\mathrm W}(G)=2$,
  $\mathrm W(G)=V\rtimes D$, $b_{\mathcal N}(G)=2$, and $H_0$ has defect three with $H_0^{\mathcal N}=R_1$ of
  defect two. The bound $b_{\mathrm{sn}}(G)\le3$ is proved in the paper and is not recomputed.

## License

The code (all `.g` and `.sh` files) is released under the MIT License (see `LICENSE`).
The output of the computations in `logs/` is released under the Creative Commons Attribution 4.0
International License (CC BY 4.0; see `logs/LICENSE`).

## How to cite

Please cite the article above. Citation metadata for this code are in `CITATION.cff`.
