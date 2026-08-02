#!/usr/bin/env sage -python
"""N609W: reconstruct the actual weighted-blowup slope form by formal sampling."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import Matrix, PolynomialRing, vector
from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def lead(x): return x[x.valuation()]
def interpolate(field, points, values, degree):
 m=Matrix(field, [[field(s)**i for i in range(degree+1)] for s in points]); return m.solve_right(vector(field,values))

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
 target,cover,pi,phi=fixture(); field,ring,u,support,matrix,xvalues=origin_data(target,cover,pi,phi); det=matrix.det(); c=lead((-support[0]/support[1])/u); formal=target.formal_group(); rows=[]
 for level in (0,1,17):
  moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); samples=[]
  for slope in [field(value) for value in range(1, 32)]:
   if slope==0 or slope==c: continue
   parameter=ring(slope)*u; point=(ring(formal.x(32)(parameter)),ring(formal.y(32)(parameter)))
   try:
    raw=sum(moving[i]*moving_basis(point,support)[i] for i in range(9))
   except ZeroDivisionError:
    continue
   normalized=-det*parameter**8*raw
   if normalized.valuation()<0: continue
   samples.append((slope, normalized[0]))
  degree=9; train=samples[:degree+1]; coeff=interpolate(field,[s for s,_ in train],[(s-c)*v for s,v in train],degree)
  R=PolynomialRing(field,'s'); s=R.gen(); numerator=sum(coeff[i]*s**i for i in range(degree+1)); heldout_ok=all((slope-c)*value==numerator(slope) for slope,value in samples[degree+1:])
  roots=numerator.roots(field); rows.append({'level':level,'sample_count':len(samples),'support_slope_constant':str(c),'numerator_degree':int(numerator.degree()),'heldout_reconstruction_match':heldout_ok,'base_field_root_count':len(roots),'base_field_multiple_root_count':sum(mult>1 for _,mult in roots),'cauchy_pole_is_not_numerator_root':numerator(c)!=0})
 gates={'three_levels_used':len(rows)==3,'all_rows_have_sufficient_slope_samples':all(row['sample_count']>=20 for row in rows),'all_degree_nine_reconstructions_replay_heldouts':all(row['heldout_reconstruction_match'] for row in rows),'all_forms_avoid_cauchy_slope_pole':all(row['cauchy_pole_is_not_numerator_root'] for row in rows)}
 o={'schema':'ecdlp.h018.poincare-direct-slope-reconstruction.n609w.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / DIRECTLY RECONSTRUCTED H018 WEIGHTED-BLOWUP SLOPE FORM / MODEL-BOUND / TOY-EVIDENCE / STRICT_TRANSFORM_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows,'normalization':'-det(V_Q)*(s*u)^8*lambda(P,Q), with v=s*u and Cauchy pole cleared by s-c'},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'The actual formal normal constants reconstruct a degree-nine slope numerator after clearing the one Cauchy slope direction and replay on all held-out base-field slopes.','next_requirement':'Use this direct slope form, not N609T, to test exceptional roots and normal derivatives, then extend to geometric roots and reciprocal charts before global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}
 a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
 if not o['preflight_pass']: raise RuntimeError('N609W direct slope reconstruction failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()
