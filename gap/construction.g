# Section 5: small instances of the construction of Theorem 5.10.
#
# Given a group E acting regularly on Omega = E, a prime p and a group
# B = V x| M_0 with a normal subgroup K of index p, put
#   X = [F_p^Omega, O^p(E)],  Y = X x| E <= C_p wr E,
# and let G be the inverse image of Y under the map B wr E -> C_p wr E
# induced by B -> B/K = C_p. In the small cases below X plays the role of the
# image iota(D) of the module D of Lemma 5.2: its dual has a vector with
# trivial stabilizer and spanning orbit.
# The script checks the conclusions of Theorem 5.10 and Proposition 5.8:
#   |G| = |K|^|E| |X| |E|,  h(G) = h(E) + h(B),  F(G) = V^Omega,
#   b_N(G) <= 2,  N = G cap B^Omega <= W(G),  l_W(G) <= l_W(E) + 1.
# Cases:
#   A. E = C_3, p = 2, B = F_7^2 x| S_3 (the shape of an odd block, but with
#      s = 2; quick);
#   D. E = C_2, p = 3, B = F_29^3 x| (C_7 x| C_3), the odd block of
#      [Part I, Proposition 6.3] (several minutes);
#   B. E = C_3, p = 2, B = AGL(2,3), a quaternion block (more than an hour;
#      run only if the variable long is set to true).
# Usage: gap -q gap/lib.g gap/construction.g
#        gap -q -c 'long:=true;;' gap/lib.g gap/construction.g

if not IsBound(long) then long := false; fi;

InverseImageConstruction := function(B, K, c, E, p)
  local n, W, emb, top, mats, X, gens, G, N;
  n := NrMovedPoints(E);
  W := WreathProduct(B, E);
  emb := List([1..n], k -> Embedding(W, k));
  top := Embedding(W, n + 1);
  mats := List(GeneratorsOfGroup(PResidual(E, p)), g -> PermutationMat(g, n, GF(p)));
  X := VectorSpace(GF(p), Concatenation(List(mats, m -> m - m^0)));
  # generators: K^Omega, the top group E, and for each basis vector f of X
  # the element of B^Omega with coordinate c^f(k) at k
  gens := Concatenation(
    Concatenation(List([1..n], k -> List(GeneratorsOfGroup(K), x -> Image(emb[k], x)))),
    List(GeneratorsOfGroup(E), e -> Image(top, e)),
    List(BasisVectors(Basis(X)), f -> Product([1..n], k -> Image(emb[k], c^IntFFE(f[k])))));
  G := Subgroup(W, gens);
  N := Intersection(G, Group(Concatenation(List([1..n],
         k -> List(GeneratorsOfGroup(B), x -> Image(emb[k], x))))));
  return rec(G := G, N := N, dimX := Dimension(X));
end;

RunCase := function(B, K, c, E, p, name)
  local r, iso, G, N, d, t;
  t := Runtime();
  r := InverseImageConstruction(B, K, c, E, p);
  iso := IsomorphismPcGroup(r.G);
  G := Image(iso); N := Image(iso, r.N);
  d := ResidualData(G);
  Print(name, ": |G| = ", Size(G),
        "; |G| = |K|^|E| |X| |E|: ", Size(G) = Size(K)^Size(E) * p^r.dimX * Size(E),
        "; h(G) = ", FittingHeight(G), " = h(E) + h(B): ", FittingHeight(G) = FittingHeight(E) + FittingHeight(B),
        "; F(G) = V^Omega: ", Size(FittingSubgroup(G)) = Size(FittingSubgroup(B))^Size(E),
        "; b_N(G) = ", d.bN,
        "; N <= W(G): ", IsSubgroup(d.W, N),
        "; |G:W(G)| = ", Size(G) / Size(d.W),
        "; l_W(G) = ", UpperWLength(G), ", l_W(E) = ", UpperWLength(E),
        "; time ", Int((Runtime() - t) / 1000), " s\n");
end;

C2 := Group((1,2));
C3 := Group((1,2,3));

# Case A
a := DiagonalMat([Z(7)^2, Z(7)^4]); b := [[0,1],[1,0]] * Z(7)^0;
G0 := SemidirectProduct(Group(a, b), GF(7)^2);
B := Image(IsomorphismPermGroup(G0)); K := DerivedSubgroup(B);
c := First(GeneratorsOfGroup(B), g -> not g in K);
RunCase(B, K, c, C3, 2, "Case A: E = C_3, p = 2, B = F_7^2 x| S_3");

# Case D
zeta := Z(29)^4;
a := DiagonalMat([zeta, zeta^2, zeta^4]); b := PermutationMat((1,2,3), 3, GF(29));
G0 := SemidirectProduct(Group(a, b), GF(29)^3);
B := Image(IsomorphismPermGroup(G0)); K := DerivedSubgroup(B);
c := First(GeneratorsOfGroup(B), g -> not g in K);
RunCase(B, K, c, C2, 3, "Case D: E = C_2, p = 3, B = F_29^3 x| (C_7 x| C_3)");

# Case B
if long then
  G0 := SemidirectProduct(GL(2,3), GF(3)^2);
  B := Image(IsomorphismPermGroup(G0)); K := DerivedSubgroup(B);
  c := First(B, g -> not g in K and Order(g) = 2);
  RunCase(B, K, c, C3, 2, "Case B: E = C_3, p = 2, B = AGL(2,3)");
fi;
Print("FINAL\n");
QUIT;
