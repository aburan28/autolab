#!/usr/bin/env sage -python
"""Independent common-regular-domain replay for N606Z."""
import argparse,json
from pathlib import Path
from sage.all import EllipticCurve,GF
import product_kummer_h018_poincare_scalar_cocycle_n606z as n606z
def main():
 p=argparse.ArgumentParser();p.add_argument("--primary",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();d=json.loads(a.primary.read_text())
 e=EllipticCurve(GF(103),[0,0,0,1,24]);q,r,s=[u for u in e.points() if not u.is_zero()][:3];sh=lambda u:84*u;vals=[]
 for x in e.points():
  if x.is_zero(): continue
  try:
   u=n606z.transport(q,r,x+s);v=n606z.transport(q+sh(r),s,x);w=n606z.transport(q,r+s,x)
   if u and v and w: vals.append(u*v/w)
  except ZeroDivisionError: pass
 c={"primary":d["preflight_pass"] is True,"samples":len(vals)>=8,"constant":len(set(vals))==1,"boundary":d["normalized_cocycle_fixed"] is False}
 o={"verified":all(c.values()),"checks":c};a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n")
 if not o["verified"]:raise RuntimeError("N606Z verifier failed")
 print(json.dumps({"verified":True,"output":str(a.output)}))
if __name__=="__main__":main()
