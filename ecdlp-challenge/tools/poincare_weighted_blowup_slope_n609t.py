#!/usr/bin/env sage -python
"""N609T: exceptional-slope leading form from the H018 weighted profile."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import PolynomialRing, vector
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def lead(x): return x[x.valuation()]

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
 target,cover,pi,phi=fixture(); field,ring,u,support,matrix,xvalues=origin_data(target,cover,pi,phi); det=matrix.det(); c=lead(support[0]*0 + (-u/u))
 # The formal support parameter is recovered from x,y: -x/y has leading c*u.
 support_parameter=-support[0]/support[1]; c=lead(support_parameter/u)
 R=PolynomialRing(field,'s'); s=R.gen(); K=R.fraction_field(); rows=[]
 for level in (0,1,17):
  moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); weighted=[-det*x for x in moving]; A=[lead(x) for x in weighted]
  poly=A[0]*s**8+A[1]*s**6-A[2]*s**5+A[3]*s**4-A[4]*s**3+A[5]*s**2-A[6]*s+A[7]
  cauchy=-A[8]*K(s)**8*((K(s)**-3+K(c)**-3)/(K(s)**-2-K(c)**-2))
  form=K(poly)+cauchy; numerator,denominator=form.numerator(),form.denominator()
  rows.append({'level':level,'support_slope_constant':int(c),'numerator_degree':int(numerator.degree()),'denominator_degree':int(denominator.degree()),'numerator_nonzero':numerator!=0,'exceptional_slope_root_count':len(numerator.roots(field)),'denominator_root_count':len(denominator.roots(field))})
 gates={'three_levels_used':len(rows)==3,'all_slope_forms_nonzero':all(row['numerator_nonzero'] for row in rows),'all_forms_have_one_genuine_cauchy_slope_pole':all(row['denominator_degree']==1 and row['denominator_root_count']==1 for row in rows)}
 o={'schema':'ecdlp.h018.poincare-weighted-blowup-slope.n609t.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / H018 WEIGHTED-BLOWUP EXCEPTIONAL SLOPE FORM MATERIALIZED / MODEL-BOUND / TOY-EVIDENCE / STRICT_TRANSFORM_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows,'substitution':'v=s*u, analyze v^8*D_t','cauchy_slope_term':'-A8*s^8*(s^-3+c^-3)/(s^-2-c^-2)'},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'The weighted pole profile produces an explicit nonzero rational slope form on the exceptional divisor, with the two Cauchy directions retained as denominator points.','next_requirement':'Compute strict-transform local equations and derivatives in both blow-up charts, then determine whether the weighted blow-up resolves the crossing before global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}
 a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
 if not o['preflight_pass']: raise RuntimeError('N609T slope form failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()
