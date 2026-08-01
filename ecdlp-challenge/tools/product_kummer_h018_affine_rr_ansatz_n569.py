#!/usr/bin/env sage -python
"""N569 exact closed-orbit affine interpolation preflight."""
import argparse,json
from pathlib import Path
from sage.all import GF,PolynomialRing,Matrix
P=103
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);a=p.parse_args();F=GF(P);R=PolynomialRing(F,"x");x=R.gen();f=x**3+x+24;u=62+102*x+93*x**2;bounds=[(0,2),(1,1),(1,2),(2,2),(2,3),(3,3),(4,4)];rows=[]
 for du,dx in bounds:
  mons=[(i,j) for i in range(du+1) for j in range(dx+1)];cols=[]
  for i,j in mons:
   q=(u**i*x**j)%f;cols.append([F(q[k]) for k in range(3)])
  M=Matrix(F,3,len(mons),lambda r,c:cols[c][r]);rows.append({"du":du,"dx":dx,"dimension":len(mons),"rank":int(M.rank()),"nullity":len(mons)-int(M.rank())})
 # H018 graph is u-(62+102x+93x^2), represented in (1,2) box; mutation values are nonzero mod f.
 graph_zero=((u-(62+102*x+93*x**2))%f)==0;f3_zero=((u-(41+10*x**2))%f)==0;identity_zero=((u-x)%f)==0
 gates={"h018_graph_control":graph_zero,"f3_mutation_rejected":not f3_zero,"identity_mutation_rejected":not identity_zero,"no_box_has_admissible_nullity_two":all(r["nullity"]!=2 for r in rows)}
 out={"schema":"ecdlp.product-kummer.h018.affine-rr-ansatz.n569.v1","claim_status":"NEGATIVE RESULT / AFFINE_NODE_INTERPOLATION_UNDERDETERMINED / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM","rows":rows,"gates":gates,"preflight_pass":all(gates.values()),"next_requirement":"Use global divisor and exceptional-node valuation data; closed-node affine interpolation alone cannot construct the H018 two-section space."};a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n");print(json.dumps({"checks":sum(gates.values()),"output":str(a.output)}))
if __name__=="__main__":main()
