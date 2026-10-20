# Section 3: residuals of affine groups G = V x| M.
#
# For each affine group in the list below the script checks:
#  - Theorem 3.3: every subnormal H of G with image A in M has
#    H^N = [V, O^p(A)] A^N. When dim V <= 4 this is checked by brute force over
#    all subnormal subgroups of G; for larger V every subspace of V is
#    subnormal and the brute force is out of reach, so the formula is checked
#    for the full inverse images VA (equation (3.4)), and b_N(G) and W(G) are
#    computed through G/F(G) as allowed by [Part I, Lemma 2.2];
#  - Theorem 3.1: b_N(G) <= 2 if and only if b_N(M) <= 2 and [V, O^p(A)] is
#    invariant under (A^N)^M for every subnormal A of M (condition (iii));
#  - Theorem 3.5: W(G) = V x| (W(M) cap the stabilizers of the [V, O^p(A)]);
#  - Proposition 3.8: b_N(M) <= b_N(G) <= b_N(M) + 1.
# Output: one line per group, then FINAL. A line with "false" is a failure.
# Usage: gap -q gap/lib.g gap/affine.g

TestAffine := function(M, p, name)
  local A, G, V, Mg, hom, reps, ok, H, Ah, U, d, repsM, bM, crit, J, A0, R0, U0, method;
  A := AffinePcGroup(M, p);
  G := A.G; V := A.V; Mg := A.M;
  hom := NaturalHomomorphismByNormalSubgroup(G, V);
  ok := true;
  if DimensionOfMatrixGroup(M) <= 4 then
    # Theorem 3.3, brute force over all subnormal subgroups of G
    method := "brute force";
    reps := SubnormalReps(G);
    for H in reps do
      Ah := PreImage(RestrictedMapping(hom, Mg), Image(hom, H));
      U := CommutatorSubgroup(V, PResidual(Ah, p));
      if NilpotentResidual(H) <> ClosureGroup(U, NilpotentResidual(Ah)) then
        ok := false;
      fi;
    od;
    d := ResidualDataBrute(G);
  else
    # Theorem 3.3 for the subgroups VA, A subnormal in M
    method := "VA only";
    reps := SubnormalReps(Mg);
    for Ah in reps do
      U := CommutatorSubgroup(V, PResidual(Ah, p));
      if NilpotentResidual(ClosureGroup(V, Ah)) <> ClosureGroup(U, NilpotentResidual(Ah)) then
        ok := false;
      fi;
    od;
    d := ResidualData(G);
  fi;
  # Theorem 3.1 (iii), Theorem 3.5, Proposition 3.8
  repsM := SubnormalReps(Mg);
  bM := 0; crit := true; J := Mg;
  for A0 in repsM do
    R0 := NilpotentResidual(A0);
    bM := Maximum(bM, Defect(Mg, R0));
    U0 := CommutatorSubgroup(V, PResidual(A0, p));
    if not IsSubgroup(Normalizer(G, U0), NormalClosure(Mg, R0)) then
      crit := false;
    fi;
    J := Intersection(J, Core(Mg, Intersection(Normalizer(Mg, R0), Normalizer(Mg, U0))));
  od;
  crit := crit and bM <= 2;
  Print(name, ": |G| = ", Size(G), ", ", method, ", classes tested = ", Length(reps),
        ", b_N(G) = ", d.bN, ", b_N(M) = ", bM,
        "; Thm 3.3: ", ok,
        "; Thm 3.1: ", crit = (d.bN <= 2),
        "; Thm 3.5: ", d.W = ClosureGroup(V, J),
        "; Prop 3.8: ", bM <= d.bN and d.bN <= bM + 1, "\n");
end;

S3wrC2 := Group((1,2,3), (1,2), (4,5,6), (4,5), (1,4)(2,5)(3,6));
C3wrC2 := Group((1,2,3), (4,5,6), (1,4)(2,5)(3,6));

TestAffine(Group([[0,1],[-1,0]] * Z(3)^0, [[1,0],[0,-1]] * Z(3)^0), 3, "D_8 on F_3^2 (Example 3.11)");
TestAffine(GL(2,3), 3, "AGL(2,3)");
TestAffine(SL(2,3), 3, "ASL(2,3)");
TestAffine(PermutationMatrices(SymmetricGroup(4), 2), 2, "S_4 on F_2^4 (permutation module)");
TestAffine(PermutationMatrices(SymmetricGroup(4), 3), 3, "S_4 on F_3^4 (permutation module)");
TestAffine(PermutationMatrices(Group((1,2,3,4), (1,3)), 3), 3, "D_8 on F_3^4 (permutation module)");
TestAffine(PermutationMatrices(RegularPermGroup(SymmetricGroup(3)), 2), 2, "S_3 on F_2^6 (regular module)");
TestAffine(PermutationMatrices(RegularPermGroup(SymmetricGroup(3)), 3), 3, "S_3 on F_3^6 (regular module)");
TestAffine(PermutationMatrices(C3wrC2, 2), 2, "C_3 wr C_2 on F_2^6 (permutation module)");
TestAffine(PermutationMatrices(S3wrC2, 2), 2, "S_3 wr C_2 on F_2^6 (permutation module)");
TestAffine(PermutationMatrices(S3wrC2, 3), 3, "S_3 wr C_2 on F_3^6 (permutation module)");
Print("FINAL\n");
QUIT;
