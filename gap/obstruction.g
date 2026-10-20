# Section 6: the obstruction in regular wreath products.
#
# Theorem 6.1, on a small instance: E = D wr F with D = S_4, F = C_2 and
# p = 3. The subgroup L = <(1,2)(3,4)> of D is nilpotent and subnormal, and
# D^N = A_4 does not normalize O^3(L) = L. For A = D_1 x L_2 the script checks
#  - that A is subnormal in E and (A^N)^E is not contained in N_E(O^3(A)),
#    so that condition (iii) of Theorem 3.10 fails for E and p = 3;
#  - that already for the permutation module V = F_3^8 of E the tested
#    residual (VA)^N has defect three in V x| E.
# Corollary 6.2: the hypotheses of Theorem 6.1 for the blocks used there:
#  (i)  B_13: no line of the module S of Q_8 over F_13 is Q_8-invariant;
#  (ii) AGL(2,3), p = 2: D^N does not normalize O^2(L) for a line L of F_3^2;
#  (iii) the block of [Part I, Proposition 6.3], p = 3: D^N does not normalize
#       O^3(L) for L = <e_1 + e_2>.
# Usage: gap -q gap/lib.g gap/obstruction.g

# ---- Theorem 6.1 ----
D := SymmetricGroup(4);
Ew := WreathProduct(D, CyclicGroup(IsPermGroup, 2));   # acts on 8 points
D1 := Image(Embedding(Ew, 1), D);
L2 := Group(Image(Embedding(Ew, 2), (1,2)(3,4)));
A := ClosureGroup(D1, L2);
Print("Theorem 6.1: A subnormal in E: ", IsSubnormal(Ew, A),
      "; (A^N)^E <= N_E(O^3(A)): ",
      IsSubgroup(Normalizer(Ew, PResidual(A, 3)), NormalClosure(Ew, NilpotentResidual(A))),
      " (expected false)\n");
PM := g -> PermutationMat(g, 8, GF(3));
G0 := SemidirectProduct(Group(List(GeneratorsOfGroup(Ew), PM)), GF(3)^8);
iso := IsomorphismPcGroup(G0);
G := Image(iso);
V := Image(iso, Image(Embedding(G0, 2)));
VA := ClosureGroup(V, Image(iso, Image(Embedding(G0, 1), Group(List(GeneratorsOfGroup(A), PM)))));
Print("Theorem 6.1: |V x| E| = ", Size(G), ", |VA| = ", Size(VA), "; defect of (VA)^N = ",
      Defect(G, NilpotentResidual(VA)), " (expected 3)\n");

# ---- Corollary 6.2 ----
cc := Z(13)^3;                                         # cc^2 = -1
Q13 := Group([[0,1],[-1,0]] * Z(13)^0, [[cc,0],[0,-cc]] * Z(13)^0);
Print("Corollary 6.2(i): Q_8-invariant lines of S over F_13: ",
      Number(Orbits(Q13, NormedRowVectors(GF(13)^2), OnLines), o -> Length(o) = 1),
      " (expected 0)\n");

B := AffinePcGroup(GL(2,3), 3);
L := Subgroup(B.G, [First(GeneratorsOfGroup(B.V), v -> not IsOne(v))]);
Print("Corollary 6.2(ii): AGL(2,3): |B:B'| = ", Index(B.G, DerivedSubgroup(B.G)),
      "; B^N normalizes O^2(L): ",
      IsSubgroup(Normalizer(B.G, PResidual(L, 2)), NilpotentResidual(B.G)), " (expected false)\n");

zeta := Z(29)^4;
a := DiagonalMat([zeta, zeta^2, zeta^4]); b := PermutationMat((1,2,3), 3, GF(29));
G0 := SemidirectProduct(Group(a, b), GF(29)^3);
L := Subgroup(G0, [Image(Embedding(G0, 2), [1,1,0] * Z(29)^0)]);
Print("Corollary 6.2(iii): |B:B'| = ", Index(G0, DerivedSubgroup(G0)),
      "; B^N normalizes O^3(L): ",
      IsSubgroup(Normalizer(G0, PResidual(L, 3)), NilpotentResidual(G0)), " (expected false)\n");
Print("FINAL\n");
QUIT;
