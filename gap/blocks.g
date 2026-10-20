# Blocks: Lemmas 4.3, 4.4, 5.4 and 5.5.
#
# Lemma 5.4: the matrices i, j, u, c over F_3 give Q = <i,j> = Q_8 = F(GL(2,3)),
#   T_0 x| C = <u,c> = S_3 complementing Q, and T_0 x| C maps isomorphically
#   onto Out(Q).
# Lemma 5.5: for the quaternion block AGL(2,3) and the odd block
#   F_29^3 x| (C_7 x| C_3) of [Part I, Proposition 6.3]: B is a W-group,
#   F(B) = V, h(B) = 4 or 3, K^N = P, some z in Q x| T_0 has C_V(z) = 0, and
#   no proper subnormal subgroup L of K satisfies LP = K.
# Lemma 4.3: for l = 13, 17 the matrices i, j with c^2 = -1 generate Q_8,
#   acting without nonzero fixed points; in V_l = Ind_Q^GL(2,3) S, a nonzero
#   vector of V_l^* supported at one induced coordinate has trivial stabilizer
#   and spans V_l^*, and these vectors form (l^2-1)/8 orbits.
# Lemma 4.4: for the odd block (s,r,t) = (3,7,29), V^* has (t-1)/r = 4
#   regular orbits of cyclic vectors supported at one coordinate.
# Usage: gap -q gap/lib.g gap/blocks.g

# ---- Lemma 5.4 ----
o := Z(3)^0;
i := [[0,1],[-1,0]] * o; j := [[1,1],[1,-1]] * o;
u := [[1,1],[0,1]] * o;  c := [[1,0],[0,-1]] * o;
M0 := GL(2,3); Q := Group(i, j); S := Group(u, c);
CQ := Centralizer(M0, Q);
Print("Lemma 5.4: Q = Q_8: ", IdGroup(Q) = [8,4],
      "; Q = F(GL(2,3)): ", Q = FittingSubgroup(M0),
      "; <u,c> = S_3: ", IdGroup(S) = [6,1],
      "; complement: ", Size(Intersection(Q, S)) = 1 and Size(ClosureGroup(Q, S)) = 48,
      "; T_0 x| C embeds in Out(Q): ", Size(Intersection(S, ClosureGroup(Q, CQ))) = 1,
      "; |Out(Q)| = 6: ", Size(AutomorphismGroup(Q)) / Size(Q / Centre(Q)) = 6, "\n");

# ---- Lemma 5.5 ----
CheckBlock := function(B, V, K, P, z, h, name)
  local d, reps, bad;
  d := ResidualData(B);
  reps := SubnormalReps(K);
  bad := Filtered(reps, L -> Size(L) < Size(K) and Size(ClosureGroup(L, P)) = Size(K));
  Print("Lemma 5.5, ", name, ": W(B) = B: ", d.W = B,
        "; F(B) = V: ", FittingSubgroup(B) = V,
        "; h(B) = ", FittingHeight(B), " (expected ", h, ")",
        "; K^N = P: ", NilpotentResidual(K) = P,
        "; C_V(z) = 1: ", Size(Centralizer(V, z)) = 1,
        "; subnormal classes of K = ", Length(reps),
        ", proper L with LP = K: ", Length(bad), "\n");
end;

A := AffinePcGroup(GL(2,3), 3);
K := DerivedSubgroup(A.G);                     # V x| SL(2,3)
P := ClosureGroup(A.V, FittingSubgroup(A.M));  # V x| Q
z := First(Centre(A.M), g -> Order(g) = 2);    # -I
CheckBlock(A.G, A.V, K, P, z, 4, "AGL(2,3)");

zeta := Z(29)^4;                               # order 7
a := DiagonalMat([zeta, zeta^2, zeta^4]); b := PermutationMat((1,2,3), 3, GF(29));
A := AffinePcGroup(Group(a, b), 29);
K := DerivedSubgroup(A.G);                     # V x| C_7
z := First(GeneratorsOfGroup(FittingSubgroup(A.M)), g -> Order(g) = 7);
CheckBlock(A.G, A.V, K, A.V, z, 3, "odd block F_29^3 x| (C_7 x| C_3)");

# ---- Lemma 4.3 ----
# Matrix of g on Ind_Q^M S, for a right transversal T of Q in M and a
# representation hom of Q on S = F_l^2.
InducedMatrix := function(T, Q, hom, g, l)
  local m, s, k, q;
  m := NullMat(2 * Length(T), 2 * Length(T), GF(l));
  for s in [1..Length(T)] do
    k := PositionCanonical(T, T[s] * g);
    q := T[s] * g * T[k]^-1;
    m{[2*s-1, 2*s]}{[2*k-1, 2*k]} := Image(hom, q);
  od;
  return ImmutableMatrix(GF(l), m);
end;

# The nonzero vectors of F_l^12 = (F_l^2)^6 supported at one coordinate.
SingleCoordinateVectors := function(l)
  local Id, vecs, s, w;
  Id := IdentityMat(12, GF(l));
  vecs := [];
  for s in [1..6] do
    for w in GF(l)^2 do
      if not IsZero(w) then Add(vecs, w[1] * Id[2*s-1] + w[2] * Id[2*s]); fi;
    od;
  od;
  return vecs;
end;

CheckLemma43 := function(l)
  local cc, ii, jj, Ql, fpf, hom, T, mats, Ml, dual, vecs, orbs, v;
  cc := Z(l)^((l-1)/4);                        # cc^2 = -1
  ii := [[0,1],[-1,0]] * Z(l)^0; jj := [[cc,0],[0,-cc]] * Z(l)^0;
  Ql := Group(ii, jj);
  fpf := ForAll(Filtered(AsList(Ql), g -> not IsOne(g)), g -> RankMat(g - One(g)) = 2);
  # induce the 2-dimensional module of Q = F(GL(2,3)) to GL(2,3)
  hom := GroupHomomorphismByImages(Q, Ql, [i, j], [ii, jj]);
  T := RightTransversal(M0, Q);
  mats := List(GeneratorsOfGroup(M0), g -> InducedMatrix(T, Q, hom, g, l));
  Ml := Group(mats);
  dual := Group(List(mats, m -> TransposedMat(m^-1)));
  vecs := SingleCoordinateVectors(l);
  orbs := Orbits(dual, vecs, OnRight);
  v := vecs[1];
  Print("Lemma 4.3, l = ", l, ": Q_8: ", IdGroup(Ql) = [8,4], "; fixed-point-free: ", fpf,
        "; |GL(2,3) on V_l| = ", Size(Ml),
        "; single-coordinate vectors: ", Length(vecs),
        ", all regular: ", ForAll(orbs, x -> Length(x) = 48),
        ", spanning: ", RankMat(Orbit(dual, v, OnRight)) = 12,
        ", orbits = ", Length(orbs), " (expected ", (l^2-1)/8, ")\n");
end;
CheckLemma43(13);
CheckLemma43(17);

# ---- Lemma 4.4 ----
Mo := Group(a, b);
dualo := Group(List(GeneratorsOfGroup(Mo), m -> TransposedMat(m^-1)));
vecs := Concatenation(List([1..3], s -> List(Filtered(AsList(GF(29)), x -> not IsZero(x)),
          x -> x * IdentityMat(3, GF(29))[s])));
orbs := Orbits(dualo, vecs, OnRight);
Print("Lemma 4.4, (s,r,t) = (3,7,29): single-coordinate vectors: ", Length(vecs),
      ", all regular: ", ForAll(orbs, x -> Length(x) = 21),
      ", all cyclic: ", ForAll(orbs, x -> RankMat(x) = 3),
      ", orbits = ", Length(orbs), " (expected ", (29-1)/7, ")\n");
Print("FINAL\n");
QUIT;
