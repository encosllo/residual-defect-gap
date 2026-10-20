# Theorem 7.10: an explicit group with b_N = 2 and b_sn = 3.
#
# G = V x| M with D = F_5^2, M = D x| <a>, a of order 3 acting by
# (x,y) -> (-y, x-y), and V = F_17^12 = (F_17^4)^3, on which (x,y) acts by
# multiplication by (zeta^y, zeta^-x, zeta^(x-y)) for zeta of order 5 in
# F_{17^4} (written as the companion matrix C of its minimal polynomial) and a
# permutes the three blocks cyclically. The script builds G as a pc group and
# checks:
#   |G| = 75 * 17^12, F(G) = V, h(G) = 3, l_W(G) = 2, W(G) = N = V x| D,
#   b_N(G) = 2 (through the tested residuals (VA)^N, A subnormal in M),
#   and, for H_0 = (R_1 + Lambda) x| Pi_1: Def_G(H_0) = 3, H_0^N = R_1 and
#   Def_G(R_1) = 2.
# The bound b_sn(G) <= 3 is proved in the paper and is not recomputed here.
# Usage: gap -q gap/lib.g gap/bryce.g

q := 17;
zeta := PrimitiveRoot(GF(q^4))^((q^4 - 1) / 5);
C := CompanionMat(MinimalPolynomial(GF(q), zeta));
Blk := function(X1, X2, X3)
  local m;
  m := NullMat(12, 12, GF(q));
  m{[1..4]}{[1..4]} := X1; m{[5..8]}{[5..8]} := X2; m{[9..12]}{[9..12]} := X3;
  return m;
end;
ActD := function(x, y) return Blk(C^y, C^(-x), C^(x - y)); end;
I4 := IdentityMat(4, GF(q));
Pa := NullMat(12, 12, GF(q));
Pa{[1..4]}{[5..8]} := I4; Pa{[5..8]}{[9..12]} := I4; Pa{[9..12]}{[1..4]} := I4;

# M as a permutation group on the points of F_5^2: translations and a
pts := Elements(GF(5)^2);
Tr := w -> PermList(List(pts, v -> Position(pts, v + w)));
amat := [[0,1],[-1,-1]] * Z(5)^0;              # (x,y) -> (-y, x-y)
ap := PermList(List(pts, v -> Position(pts, v * amat)));
e1 := [1,0] * Z(5)^0; e2 := [0,1] * Z(5)^0;
M := Group(Tr(e1), Tr(e2), ap);
rho := GroupHomomorphismByImages(M, GL(12, GF(q)), GeneratorsOfGroup(M),
                                 [ActD(1,0), ActD(0,1), Pa]);
Print("|M| = ", Size(M), "; the action is a homomorphism: ", rho <> fail, "\n");

pciso := IsomorphismPcGroup(M); Mpc := Image(pciso);
pcg := Pcgs(Mpc);
modmats := List(pcg, g -> Image(rho, PreImagesRepresentative(pciso, g)));
G := SplitExtension(Mpc, GModuleByMats(modmats, GF(q)));
pp := Pcgs(G); n0 := Length(pcg);
Vgens := pp{[n0+1..n0+12]};
V := Subgroup(G, Vgens);
Comp := Subgroup(G, pp{[1..n0]});
emb := GroupHomomorphismByImages(Mpc, Comp, AsList(pcg), pp{[1..n0]});
Dsub := Image(emb, FittingSubgroup(Mpc));
N := ClosureGroup(V, Dsub);
Print("|G| = 75 * 17^12: ", Size(G) = 75 * 17^12,
      "; F(G) = V: ", FittingSubgroup(G) = V,
      "; h(G) = ", FittingHeight(G), "\n");

# tested residuals: (VA)^N for A subnormal in M (Part I, Lemma 2.2)
reps := SubnormalReps(Mpc);
res := List(reps, A -> NilpotentResidual(ClosureGroup(V, Image(emb, A))));
W := G;
for R in res do W := Intersection(W, Core(G, Normalizer(G, R))); od;
Print("subnormal classes of M = ", Length(reps),
      "; defects of the tested residuals = ", List(res, R -> Defect(G, R)),
      "; W(G) = N: ", W = N,
      "; l_W(G) = ", UpperWLength(G), "\n");

# H_0 = (R_1 + Lambda) x| Pi_1, with Pi_1 = {(x,0)}, R_1 = 0 + V_2 + V_3 and
# Lambda the line spanned by the first basis vector of V_1
Pi1 := Subgroup(G, [Image(emb, Image(pciso, Tr(e1)))]);
Print("Pi_1 centralizes V_1: ", ForAll([1..4], k -> Vgens[k]^GeneratorsOfGroup(Pi1)[1] = Vgens[k]), "\n");
R1 := Subgroup(G, Vgens{[5..12]});
H0 := ClosureGroup(ClosureGroup(R1, Vgens[1]), Pi1);
Print("|H_0| = 17^9 * 5: ", Size(H0) = 17^9 * 5,
      "; Def_G(H_0) = ", Defect(G, H0),
      "; H_0^N = R_1: ", NilpotentResidual(H0) = R1,
      "; Def_G(R_1) = ", Defect(G, R1), "\n");
Print("FINAL\n");
QUIT;
