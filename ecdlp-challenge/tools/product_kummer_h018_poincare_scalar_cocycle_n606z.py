#!/usr/bin/env sage -python
"""N606Z: explicit Miller scalar transport and cocycle test for N606Y."""
import argparse, json
from pathlib import Path
from sage.all import EllipticCurve, GF

def line_value(p, q, x):
    if p + q == p.curve()(0): return x[0] - p[0]
    if p == q:
        slope = (3*p[0]**2 + p.curve().a4()) / (2*p[1])
    else: slope = (q[1]-p[1])/(q[0]-p[0])
    r = p + q
    return (x[1]-p[1]-slope*(x[0]-p[0]))/(x[0]-r[0])

def miller_multiple(p, count, x):
    current, value = p, x[0].parent()(1)
    for bit in bin(count)[3:]:
        value = value**2 * line_value(current, current, x); current = current + current
        if bit == "1": value *= line_value(current, p, x); current += p
    return value

def transport(q, r, x):
    zq = (-4)*q; a = zq-r
    return miller_multiple(-r, 8, x) * line_value((-8)*r, a, x)

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    e=EllipticCurve(GF(103),[0,0,0,1,24]); pts=[u for u in e.points() if not u.is_zero()]
    q,r,s=pts[0],pts[1],pts[2]; shift=lambda u:84*u
    ratios=[]
    for x in pts:
        try:
            first=transport(q,r,x+s); second=transport(q+shift(r),s,x); right=transport(q,r+s,x)
            if first and second and right: ratios.append(first*second/right)
        except ZeroDivisionError: pass
    gates={"regular_samples":len(ratios)>=8,"cocycle_ratio_constant":len(set(ratios))==1,"base_shift":(-4)*(q+shift(r))==(-4)*q-9*r}
    out={"schema":"ecdlp.product-kummer.h018.poincare-scalar-cocycle.n606z.v1","claim_status":"OBSERVATION / EXPLICIT_MILLER_SCALAR_TRANSPORT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM","records":{"common_regular_samples":len(ratios),"cocycle_constant":str(ratios[0]) if ratios else None},"gates":gates,"preflight_pass":all(gates.values()),"normalized_cocycle_fixed":False,"explicit_bundle_isomorphism_constructed":False,"next_requirement":"Fix the remaining scalar normalization against a Poincare rigidification and compare the resulting frame with N606U."}
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="ascii")
    if not out["preflight_pass"]: raise RuntimeError("N606Z failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(a.output)},sort_keys=True))
if __name__=="__main__": main()
