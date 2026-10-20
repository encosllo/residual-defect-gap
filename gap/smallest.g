# Example 3.11: the smallest group with residual defect two.
#
# Claims checked:
#  (a) every group of order less than 72 has b_N <= 1;
#  (b) S_3 wr C_2 = SmallGroup(72,40) is the only group of order 72 with
#      b_N = 2, and its ordinary defect is b_sn = 2;
#  (c) in the affine model G = V x| M with V = F_3^2 and M = D_8 the group of
#      signed permutation matrices, b_N(G) = 2 and W(G) = V x| Z(M), so that
#      |W(G)| = 18 and |G:W(G)| = 4.
# Nilpotent groups are skipped: every tested residual is trivial, so b_N <= 1.
# Usage: gap -q gap/lib.g gap/smallest.g

found := [];
for n in [1..72] do
  for i in [1..NrSmallGroups(n)] do
    G := SmallGroup(n, i);
    if not IsNilpotentGroup(G) then
      d := ResidualDataBrute(G);
      if d.bN >= 2 then Add(found, [n, i, d.bN, d.bsn]); fi;
    fi;
  od;
od;
Print("groups of order at most 72 with b_N >= 2 [n, i, b_N, b_sn]: ", found, "\n");

W := WreathProduct(SymmetricGroup(3), CyclicGroup(IsPermGroup, 2));
Print("IdGroup(S_3 wr C_2) = ", IdGroup(W), "\n");

x := [[0,1],[-1,0]] * Z(3)^0;
y := [[1,0],[0,-1]] * Z(3)^0;
A := AffinePcGroup(Group(x, y), 3);
d := ResidualDataBrute(A.G);
Print("affine model V x| D_8: |G| = ", Size(A.G), ", IdGroup = ", IdGroup(A.G),
      ", b_N = ", d.bN, ", |W(G)| = ", Size(d.W),
      ", W(G) = V x| Z(M): ", d.W = ClosureGroup(A.V, Centre(A.M)), "\n");
Print("FINAL\n");
QUIT;
