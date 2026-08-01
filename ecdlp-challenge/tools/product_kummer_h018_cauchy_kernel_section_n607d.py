#!/usr/bin/env sage -python
"""N607D moving Cauchy-kernel fibre-basis preflight."""
import argparse,json
from pathlib import Path
from sage.all import EllipticCurve,GF,Matrix
def vals(p,r,a):
 x,y=p[0],p[1]; xr,yr=r[0],r[1]
 k=(y+yr)/(x-xr) if x!=xr else None
 return [1,x,y,x*x,x*y,x**3,x*x*y,x**4,k]
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);a=p.parse_args();E=EllipticCurve(GF(103),[0,0,0,1,24]);qs=[u for u in E.points() if not u.is_zero()][:3];rows=[]
 for q in qs:
  r=(-4)*q; candidates=[u for u in E.points() if not u.is_zero() and u not in (r,-r)][:9];M=Matrix(GF(103),[vals(u,r,E.a4()) for u in candidates]);limit=-(3*r[0]**2+E.a4())/(2*r[1]);rows.append({"rank":int(M.rank()),"minus_r_limit":str(limit),"r_y_nonzero":r[1]!=0})
 gates={"all_fibre_ranks_nine":all(x["rank"]==9 for x in rows),"minus_r_regular_limit":all(x["r_y_nonzero"] for x in rows),"sign_flipped_has_wrong_pole":all(((-r)[1]-r[1])!=0 for r in [(-4)*q for q in qs])}
 o={"schema":"ecdlp.product-kummer.h018.cauchy-kernel-section.n607d.v1","claim_status":"OBSERVATION / MOVING_CAUCHY_KERNEL_FIBRE_BASIS / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM","rows":rows,"gates":gates,"preflight_pass":all(gates.values()),"global_h018_sections_constructed":False,"next_requirement":"Impose N606Y/N606Z transport compatibility on coefficient functions to determine whether two global H018 sections exist in this basis."};a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps({"checks":sum(gates.values()),"output":str(a.output)}))
if __name__=="__main__":main()
