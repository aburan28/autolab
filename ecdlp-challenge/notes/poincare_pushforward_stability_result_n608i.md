# Result: N608I H018 Abstract Pushforward Identification

## Claim

On the registered `F_103` H018 fixture, the two rank-nine bundles

```text
F = p2_* (p1^* O(9O) tensor (id,z)^* P tensor p2^* O(11O))
G = pi_* O_Eprime(2Oprime)
```

are isomorphic as abstract vector bundles, provided the standard
Fourier--Mukai and elliptic-bundle classification theorems listed below.  This
upgrades N606U's unchecked stability step; it does not construct an isomorphism
matrix.

## Exact Fixture Evidence

- `z=-3-pi` has degree `97` and is separable in characteristic `103`.
- The normalized Fourier--Mukai bundle `W9` has rank `9`, degree `-1`, and
  determinant `O(-O)`.
- Therefore `z^*W9` has degree `-97`; tensoring by `O(11O)` adds `9*11=99`.
  Thus `F` has rank `9`, degree `2`, and determinant `O(2O)`.
- Separable pullback makes `z^*W9` semistable; `gcd(9,97)=1` makes it stable,
  and tensoring by `O(11O)` preserves stability.  Hence `F` is stable.
- The N606U cover `pi:Eprime -> E` has cyclic degree `9`.  Its order-nine
  kernel generator is explicit over `F_(103^6)`.
- Pulling `G` back along `pi` gives the nine translates of `O(2Oprime)`.  A
  nonzero kernel point cannot stabilize this line bundle because that would
  require it to be two-torsion, while the kernel has odd order `9`.
- Hence the pullback is a direct sum of degree-two line bundles and is
  semistable.  Semistability descends through the finite separable cover, so
  `G` is semistable.  Since `gcd(9,2)=1`, it is stable.
- The norm gives `Nm(O(2Oprime))=O(2O)`.  The determinant of `pi_*O` is
  trivial: over the cyclic deck group it is the product of all nine characters,
  whose exponent sum is `0+...+8=36=0 mod 9`.  Therefore
  `det(G)=O(2O)`.

The producer and independently recomputed arithmetic both pass:

- Contract: `notes/poincare_pushforward_stability_contract_n608i.md`
- Producer: `tools/poincare_pushforward_stability_n608i.py`
- Receipt: `notes/poincare_pushforward_stability_n608i.json`
- Independent replay: `tools/verify_poincare_pushforward_stability_n608i.py`

## Standard-Theorem Dependencies

1. Fourier--Mukai transforms on an elliptic curve take `O(9O)` to the stated
   stable rank-nine, degree-minus-one bundle.
2. Separable isogeny pullback preserves semistability; coprime rank and degree
   promote semistability to stability.
3. Stable elliptic bundles with coprime rank and degree are uniquely determined
   by rank, degree, and determinant.

The classification statement is the classical elliptic-bundle theorem of
[Atiyah](https://londmathsoc.onlinelibrary.wiley.com/doi/abs/10.1112/plms/s3-7.1.414).
A Fourier--Mukai presentation of the coprime rank-degree/determinant result is
also given by [Ploog](https://www.researchgate.net/publication/242357988_Fourier-Mukai_Transforms_and_Stable_Bundles_on_Elliptic_Curves).

## Consequence

The N606U target is no longer merely a numerical resemblance: under those
standard theorems it is the unique abstract bundle with the recorded
invariants.  The unresolved work is fully constructive:

```text
Construct M on declared base charts with
M: F -> pi_* O_Eprime(2Oprime),
verify normalized-Poincare compatibility, and then evaluate the two global
surface sections through M.
```

## Boundary

This is a `RESTRICTED THEOREM / STANDARD-THEOREM-BOUND / TOY-EVIDENCE`
existence result. It does not supply `M`, a section evaluator, a factor base,
source recovery, relation rank, target descent, a charged comparison with rho,
or an ECDLP improvement.
