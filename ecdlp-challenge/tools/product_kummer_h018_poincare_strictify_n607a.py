#!/usr/bin/env sage -python
"""N607A: extension-field cyclic strictification of N606Z transports."""
import argparse,json
from math import gcd
from pathlib import Path
from sage.all import EllipticCurve,GF,inverse_mod
import product_kummer_h018_poincare_scalar_cocycle_n606z as z
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);a=p.parse_args();E=EllipticCurve(GF(103),[0,0,0,1,24]);K=GF(103**2,name="a");EK=E.base_extend(K);g0=next(u for u in E.points() if not u.is_zero());g=EK(K(g0[0]),K(g0[1]));qs=[i*(84*g) for i in range(109)]; sample=None
 for X in EK.points():
  if X.is_zero():continue
  try:
   values=[z.transport(qs[i],g,X+(108-i)*g) for i in range(109)]
   if all(values):sample=(X,values);break
  except ZeroDivisionError:pass
 if sample is None:raise RuntimeError("no common regular extension point")
 X,values=sample;full=K.prod(values);units=K.order()-1;scale=full**(-inverse_mod(109,units));gates={"common_regular_extension_point":True,"coprime_cycle_root":gcd(109,units)==1,"full_cycle_strictified":scale**109*full==1}
 out={"schema":"ecdlp.product-kummer.h018.poincare-strictify.n607a.v1","claim_status":"OBSERVATION / EXTENSION_FIELD_CYCLIC_STRICTIFICATION / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM","records":{"field_degree":2,"field_unit_order":int(units),"cycle_order":109,"regular_point":str(X)},"gates":gates,"preflight_pass":all(gates.values()),"geometric_normalized_poincare_rigidification_constructed":False,"explicit_bundle_isomorphism_constructed":False,"next_requirement":"Prove this finite cyclic gauge descends from a geometric normalized-Poincare rigidification and compare it to N606U."};a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n");print(json.dumps({"checks":sum(gates.values()),"output":str(a.output)},sort_keys=True))
if __name__=="__main__":main()
