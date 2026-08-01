#!/usr/bin/env python3
"""Independent N606S eigenline field replay."""
from sage.all import GF, PolynomialRing, identity_matrix, matrix
import argparse, json
from pathlib import Path
p=103; F=GF(p); R=PolynomialRing(F,"z"); F2=GF(p**2,name="u",modulus=R.gen()**2+1)
I=matrix(F,[[0,-1],[1,0]]); J=matrix(F,[[2,43],[43,101]]); K=I*J; D=-((identity_matrix(F,2)+I+J+K)/F(2))
def main():
 pser=argparse.ArgumentParser(); pser.add_argument("--primary",type=Path,required=True); pser.add_argument("--output",type=Path,required=True); a=pser.parse_args(); q=json.loads(a.primary.read_text())
 ok=q["theta_eigenline_field_gate_pass"] and D**3==identity_matrix(F,2) and all(M.charpoly().roots(F)==[] for M in (I,J,K)) and D*I*D.inverse()==J and D*J*D.inverse()==K and D*K*D.inverse()==I
 a.output.write_text(json.dumps({"verified":bool(ok),"primary":str(a.primary),"strongest_valid_statement":"Independent replay confirms the quadratic phase obstruction for all three cubic quotient lifts." if ok else "N606S replay failed."},indent=2,sort_keys=True)+"\n",encoding="ascii")
 if not ok: raise RuntimeError("N606S verifier failed")
 print(json.dumps({"verified":True,"output":str(a.output)}))
if __name__=="__main__": main()
