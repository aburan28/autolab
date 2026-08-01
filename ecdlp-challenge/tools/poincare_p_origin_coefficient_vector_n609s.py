#!/usr/bin/env sage -python
"""N609S: Q-valuations of determinant-weighted raw P-origin coefficients."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import vector
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture
from poincare_origin_boundary_n609c import constant_coefficient

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def val(x): return None if x.is_zero() else int(x.valuation())

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
 target,cover,pi,phi=fixture(); field,ring,_u,_s,matrix,xvalues=origin_data(target,cover,pi,phi); det=matrix.det(); rows=[]
 for level in (0,1,17):
  moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); weighted=[-det*x for x in moving]
  rows.append({'level':level,'weighted_valuations':[val(x) for x in weighted],'weighted_constants':[int(constant_coefficient(x,field)) if val(x) is not None and val(x)>=0 else None for x in weighted]})
 # indices 1..7 are polynomial P-pole coefficients. Index 8 is Cauchy and requires its own chart.
 expected=[-8,-6,-5,-4,-3,-2,-1,0,-7]
 gates={'three_levels_used':len(rows)==3,'all_levels_have_exact_weighted_crossing_profile':all(row['weighted_valuations']==expected for row in rows),'x4_full_coefficient_retains_n609r_constant':all(row['weighted_constants'][7]==31 for row in rows)}
 o={'schema':'ecdlp.h018.poincare-p-origin-coefficient-vector.n609s.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'NEGATIVE RESULT / NAIVE PRODUCT-CHART H018 CROSSING EXTENSION OBSTRUCTED BY WEIGHTED POLE PROFILE / MODEL-BOUND / TOY-EVIDENCE / WEIGHTED_BLOWUP_CHART_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows,'polynomial_indices':[1,2,3,4,5,6,7],'cauchy_index':8,'expected_weighted_valuation_profile':expected},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'The exact coefficient profile has Q-poles matching the P-pole weights, so a naive product-chart extension is obstructed. The profile instead supplies a concrete weighted-blow-up candidate with slope coordinate v/u.','next_requirement':'Construct the weighted blow-up chart at P=Q=O, substitute v=s*u (and its reciprocal chart), and test whether the determinant has a regular nonzero strict-transform equation before global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}
 a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
 if not o['preflight_pass']: raise RuntimeError('N609S coefficient vector failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()
