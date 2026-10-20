# Theorem 3.10: the uniform criterion through the regular module.
#
# For a group M and a prime p, condition (iii) of Theorem 3.10 reads:
# b_N(M) <= 2 and (A^N)^M <= N_M(O^p(A)) for every subnormal A of M.
# Theorem 3.10 states that it holds if and only if b_N(F_p[M] x| M) <= 2.
# For each pair (M, p) below the script computes both sides.
# Usage: gap -q gap/lib.g gap/uniform.g

UniformCondition := function(M, p)
  local reps, bM;
  reps := SubnormalReps(M);
  bM := Maximum(List(reps, A -> Defect(M, NilpotentResidual(A))));
  return bM <= 2 and ForAll(reps, A ->
    IsSubgroup(Normalizer(M, PResidual(A, p)), NormalClosure(M, NilpotentResidual(A))));
end;

TestUniform := function(M, p, name)
  local R, A, d, u;
  R := RegularPermGroup(M);
  A := AffinePcGroup(PermutationMatrices(R, p), p);
  d := ResidualData(A.G);
  u := UniformCondition(M, p);
  Print(name, ", p = ", p, ": |F_p[M] x| M| = ", Size(A.G),
        ", b_N = ", d.bN, "; condition (iii): ", u,
        "; Thm 3.10: ", u = (d.bN <= 2), "\n");
end;

S3 := SymmetricGroup(3);
D8 := Group((1,2,3,4), (1,3));
A4 := AlternatingGroup(4);
C3wrC2 := Group((1,2,3), (4,5,6), (1,4)(2,5)(3,6));
TestUniform(S3, 2, "S_3");
TestUniform(S3, 3, "S_3");
TestUniform(D8, 3, "D_8");
TestUniform(A4, 2, "A_4");
TestUniform(A4, 3, "A_4");
TestUniform(C3wrC2, 2, "C_3 wr C_2");
TestUniform(C3wrC2, 3, "C_3 wr C_2");
Print("FINAL\n");
QUIT;
