#!/usr/bin/env sage -python
"""N609V: formal normal-direction test at rational H018 blow-up roots."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import PolynomialRing, vector
from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def lead(x): return x[x.valuation()]

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
 target,cover,pi,phi=fixture(); field,ring,u,support,matrix,xvalues=origin_data(target,cover,pi,phi); det=matrix.det(); support_parameter=-support[0]/support[1]; c=lead(support_parameter/u); R=PolynomialRing(field,'s'); s=R.gen(); K=R.fraction_field(); rows=[]
 for level in (0,1,17):
  moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); A=[lead(-det*x) for x in moving]
  poly=A[0]*s**8+A[1]*s**6-A[2]*s**5+A[3]*s**4-A[4]*s**3+A[5]*s**2-A[6]*s+A[7]
  form=K(poly)-A[8]*K(s)**8*((K(s)**-3+K(c)**-3)/(K(s)**-2-K(c)**-2)); num,den=form.numerator(),form.denominator()
  for root,_mult in num.roots(field):
   if den(root)==0: continue
   parameter=ring(root)*u; formal=target.formal_group(); point=(ring(formal.x(32)(parameter)),ring(formal.y(32)(parameter)))
   raw=sum(moving[i]*moving_basis(point,support)[i] for i in range(9)); normal=-det*parameter**8*raw
   rows.append({'level':level,'slope':str(root),'normal_valuation':int(normal.valuation()),'normal_leading_coefficient':str(normal[normal.valuation()])})
 gates={'all_expected_rational_exceptional_roots_tested':len(rows)==15,'purported_slope_roots_do_not_zero_direct_normal_substitution':all(row['normal_valuation']==0 for row in rows)}
 o={'schema':'ecdlp.h018.poincare-weighted-blowup-normal-probe.n609v.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'NEGATIVE RESULT / DERIVED H018 WEIGHTED-SLOPE FORM DOES NOT MATCH DIRECT NORMAL SUBSTITUTION / MODEL-BOUND / TOY-EVIDENCE / TWO_PARAMETER_ASYMPTOTICS_REDERIVATION_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'At all fifteen purported exceptional roots of the derived slope form, direct formal substitution gives a nonzero normalized constant term. Therefore the slope formula is not yet a valid strict-transform equation and must not support a blow-up smoothness or global-Cartier claim.','next_requirement':'Recompute the two-parameter Cauchy/support asymptotics directly, with a verified normalization convention, before attempting another blow-up chart, global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claim.'}
 a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
 if not o['preflight_pass']: raise RuntimeError('N609V normal probe failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()
