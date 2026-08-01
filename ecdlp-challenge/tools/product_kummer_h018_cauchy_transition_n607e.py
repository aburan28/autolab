#!/usr/bin/env sage -python
"""N607E: transition matrices for the N607D moving fibre basis."""
import argparse,json,sys
from pathlib import Path
from sage.all import EllipticCurve,GF,Matrix,identity_matrix
sys.path.insert(0,str(Path(__file__).resolve().parent))
import product_kummer_h018_poincare_scalar_cocycle_n606z as z
def basis(p,r):
 if r.is_zero():
  x,y=p[0],p[1];return [1,x,y,x*x,x*y,x**3,x*x*y,x**4,x**3*y]
 x,y=p[0],p[1];xr,yr=r[0],r[1];return [1,x,y,x*x,x*y,x**3,x*x*y,x**4,(y+yr)/(x-xr)]
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);a=p.parse_args();E=EllipticCurve(GF(103),[0,0,0,1,24]);K=GF(103**2,name="a");EK=E.base_extend(K);g0=next(u for u in E.points() if not u.is_zero());g=EK(K(g0[0]),K(g0[1]));qs=[i*(84*g) for i in range(109)];points=[u for u in EK.points() if not u.is_zero()];M=identity_matrix(K,9);transitions=0
 def row_pair(q,r,rp,x):
  f=z.transport(q,g,x)
  if f==0 or (not r.is_zero() and x in (r,-r)) or (not rp.is_zero() and x in (rp,-rp)): raise ZeroDivisionError
  return basis(x,rp),[f*v for v in basis(x+g,r)]
 for i,q in enumerate(qs):
  qp=qs[(i+1)%109];r=(-4)*q;rp=(-4)*qp;chosen=[];old=[];new=[]
  # One rational point plus four conjugate pairs gives a Frobenius-stable 9-point frame.
  for x0 in E.points():
   if x0.is_zero(): continue
   try:
    x=EK(K(x0[0]),K(x0[1]));row_new,row_old=row_pair(q,r,rp,x);chosen.append(x);new.append(row_new);old.append(row_old);break
   except ZeroDivisionError: continue
  seen=set()
  for x in points:
   if len(chosen)>=9: break
   phi=EK(x[0]**103,x[1]**103)
   if x==phi or str(x) in seen or str(phi) in seen: continue
   try:
    rn,ro=row_pair(q,r,rp,x);pn,po=row_pair(q,r,rp,phi)
    chosen.extend([x,phi]);new.extend([rn,pn]);old.extend([ro,po]);seen.update([str(x),str(phi)])
   except ZeroDivisionError: continue
  if len(chosen)!=9:raise RuntimeError("no regular interpolation set")
  C=Matrix(K,new).solve_right(Matrix(K,old));M=C*M;transitions+=1
 scalar=all(M[i,j]==(M[0,0] if i==j else 0) for i in range(9) for j in range(9))
 eigenspaces=[]
 for root,multiplicity in M.charpoly().roots(K): eigenspaces.append({"eigenvalue":str(root),"algebraic_multiplicity":int(multiplicity),"geometric_multiplicity":9-int((M-root*identity_matrix(K,9)).rank())})
 roots=[root for root,_ in M.charpoly().roots(K)]
 conjugate_pair=len(roots)==2 and roots[0]**103==roots[1] and roots[1]**103==roots[0]
 base_entries=all(value**103==value for value in M.list())
 out={"schema":"ecdlp.product-kummer.h018.cauchy-transition.n607e.v1","claim_status":"OBSERVATION / FROBENIUS_STABLE_MOVING_FIBRE_MONODROMY / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM","records":{"field_degree":2,"transition_count":transitions,"monodromy_entries_in_base_field":base_entries,"monodromy_rank":int(M.rank()),"monodromy_is_identity":M==identity_matrix(K,9),"monodromy_is_scalar":scalar,"monodromy_scalar":str(M[0,0]) if scalar else None,"eigenspaces":eigenspaces,"linear_eigenvalues_are_frobenius_conjugate":conjugate_pair},"preflight_pass":transitions==109 and M.rank()==9 and base_entries,"theta11_twist_included":False,"global_h018_sections_constructed":False,"next_requirement":"Construct a geometric theta transition and compare its algebraic global sections to this finite-rational-orbit spectrum; one-dimensional orbit eigenspaces alone do not determine H0."};a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n");print(json.dumps({"transitions":transitions,"output":str(a.output)}))
if __name__=="__main__":main()
