#!/usr/bin/env sage -python
"""N607E: transition matrices for the N607D moving fibre basis."""
import argparse,json,sys
from pathlib import Path
from sage.all import EllipticCurve,GF,Matrix,identity_matrix
sys.path.insert(0,str(Path(__file__).resolve().parent))
import product_kummer_h018_poincare_scalar_cocycle_n606z as z
def basis(p,r):
 x,y=p[0],p[1];xr,yr=r[0],r[1];return [1,x,y,x*x,x*y,x**3,x*x*y,x**4,(y+yr)/(x-xr)]
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);a=p.parse_args();E=EllipticCurve(GF(103),[0,0,0,1,24]);K=GF(103**2,name="a");EK=E.base_extend(K);g0=next(u for u in E.points() if not u.is_zero());g=EK(K(g0[0]),K(g0[1]));qs=[i*(84*g) for i in range(109)];points=[u for u in EK.points() if not u.is_zero()];M=identity_matrix(K,9);transitions=0
 for i,q in enumerate(qs):
  qp=qs[(i+1)%109];r=(-4)*q;rp=(-4)*qp;chosen=[];old=[];new=[]
  for x in points:
   try:
    f=z.transport(q,g,x)
    if f==0 or x in (r,-r,rp,-rp):continue
    row_new=basis(x,rp);row_old=[f*v for v in basis(x+g,r)]
    if None not in row_new and None not in row_old:chosen.append(x);new.append(row_new);old.append(row_old)
    if len(chosen)==9 and Matrix(K,new).rank()==9:break
   except ZeroDivisionError:pass
  if len(chosen)!=9:raise RuntimeError("no regular interpolation set")
  C=Matrix(K,new).solve_right(Matrix(K,old));M=C*M;transitions+=1
 out={"schema":"ecdlp.product-kummer.h018.cauchy-transition.n607e.v1","claim_status":"OBSERVATION / MOVING_FIBRE_TRANSITION_MONODROMY / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM","records":{"field_degree":2,"transition_count":transitions,"monodromy_rank":int(M.rank()),"monodromy_is_identity":M==identity_matrix(K,9)},"preflight_pass":transitions==109 and M.rank()==9,"theta11_twist_included":False,"global_h018_sections_constructed":False,"next_requirement":"Incorporate the degree-eleven second-factor theta transition before interpreting monodromy fixed vectors as global sections."};a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n");print(json.dumps({"transitions":transitions,"output":str(a.output)}))
if __name__=="__main__":main()
