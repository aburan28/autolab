#!/usr/bin/env sage -python
"""N609U: exceptional-divisor smoothness screen for the H018 blow-up form."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import PolynomialRing, gcd, vector
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
  form=K(poly)-A[8]*K(s)**8*((K(s)**-3+K(c)**-3)/(K(s)**-2-K(c)**-2))
  num,den=form.numerator(),form.denominator()
  common=gcd(num,num.derivative()); pole_common=gcd(num,den)
  rows.append({'level':level,'numerator_degree':int(num.degree()),'squarefree':common.degree()==0,'numerator_denominator_coprime':pole_common.degree()==0,'rational_root_count':len(num.roots(field)),'rational_multiple_root_count':sum(m>1 for _,m in num.roots(field))})
 gates={'three_levels_used':len(rows)==3,'all_exceptional_numerators_squarefree':all(row['squarefree'] for row in rows),'no_exceptional_root_hits_cauchy_pole':all(row['numerator_denominator_coprime'] for row in rows),'all_rational_exceptional_roots_simple':all(row['rational_multiple_root_count']==0 for row in rows)}
 o={'schema':'ecdlp.h018.poincare-weighted-blowup-exceptional-smoothness.n609u.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / RATIONAL EXCEPTIONAL-DIVISOR SMOOTHNESS SCREEN FOR H018 WEIGHTED BLOWUP / MODEL-BOUND / TOY-EVIDENCE / FULL_STRICT_TRANSFORM_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows,'support_slope_constant':int(c)},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'The exceptional slope numerators are tested for repeated roots and collision with the genuine Cauchy slope pole. This screens the exceptional divisor only, not the full blow-up charts.','next_requirement':'Derive the first strict-transform coefficient off the exceptional divisor and test both affine blow-up charts, then extend to non-rational points before global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}
 a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
 if not o['preflight_pass']: raise RuntimeError('N609U exceptional smoothness failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()
