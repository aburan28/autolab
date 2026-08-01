# Result: N606R Principal-Theta Quotient Descent Gate

## Status

RESTRICTED THEOREM / NEGATIVE BASE_FIELD_PRINCIPAL_THETA_QUOTIENT_SHORTCUT /
CUBIC_ORBIT_DESCENT_CONSTRUCTION_OPEN / MODEL-BOUND / TOY-EVIDENCE /
LITERATURE_BOUND / NO_ECDLP_CLAIM

## Claim

Test whether a maximal isotropic order-two subgroup of the type-(1,2) H018
kernel is defined over `F_103`.  Under the stated polarization-descent
assumption, such a subgroup would give a principally polarized quotient and a
single base-field principal theta divisor whose pullback is an H018 section.

The descent criterion is literature-derived: Lubicz and Robert describe the
isotropic-kernel descent condition and the principal-polarization conclusion
for a maximal isotropic subgroup in their 2015 isogeny construction.  This
receipt verifies only its finite-kernel field-of-definition prerequisite.

## Exact Result

The H018 kernel from N602A has three nonzero order-two lines:

```text
G_g, G_h, G_{g+h}.
```

Coordinatewise Frobenius multiplication by `pi` in `F_4` permutes them as

```text
G_g -> G_h -> G_{g+h} -> G_g.
```

Therefore no maximal isotropic line is `F_103`-stable.  The zero subgroup and
the full four-element kernel are stable, but neither is a maximal isotropic
order-two subgroup.  Frobenius cubed fixes every line, so each candidate
principal quotient is defined over `F_103^3` in this kernel model.  The
identity-Frobenius positive control stabilizes all three lines.

The primary producer passes `6/6` checks and the independent finite-kernel
replay is verified.

## Interpretation

This rules out only the shortcut of constructing a single base-field principal
theta quotient from H018.  It does not rule out H018 itself, its class-level
two-section budget, an explicit quotient over `F_103^3`, semilinear descent of
three conjugate theta divisors, a rigidified `F_103` basis, a refined-cover
evaluator, or any later factor-base, rank, descent, or ECDLP result.

## Next Concrete Action

Construct the quotient by `G_g` over `F_103^3`, materialize its principal theta
divisor and the two Frobenius conjugates, then compute their pullbacks to H018.
Only an explicit semilinear descent calculation can turn that orbit into named
`F_103` scalar sections `s0,s1`.
