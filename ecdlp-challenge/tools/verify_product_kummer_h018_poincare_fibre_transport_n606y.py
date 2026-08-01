#!/usr/bin/env sage -python
"""Independent replay for N606Y."""
import argparse, json
from pathlib import Path
from sage.all import EllipticCurve, GF

def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--primary",type=Path,required=True); parser.add_argument("--output",type=Path,required=True); args=parser.parse_args()
    primary=json.loads(args.primary.read_text(encoding="ascii")); curve=EllipticCurve(GF(103),[0,0,0,1,24]); pts=[p for p in curve.points() if not p.is_zero()]; q,r=pts[:2]
    checks={"primary_pass":primary["preflight_pass"] is True,"shift_identity":(-4)*(q+84*r)==(-4)*q-9*r,"dimension":len(curve.riemann_roch_basis(8*curve.divisor(curve(0))+curve.divisor((-4)*q)))==9,"boundary":primary["scalar_rigidification_constructed"] is False and primary["explicit_bundle_isomorphism_constructed"] is False}
    out={"verified":all(checks.values()),"checks":checks,"strongest_valid_statement":"Independent replay confirms only the discrete Poincare divisor-class transport, not a scalar cocycle or surface evaluator."}; args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="ascii")
    if not out["verified"]: raise RuntimeError("N606Y verifier failed")
    print(json.dumps({"verified":True,"output":str(args.output)},sort_keys=True))
if __name__=="__main__": main()
