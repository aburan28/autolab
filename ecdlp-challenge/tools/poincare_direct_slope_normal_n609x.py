#!/usr/bin/env sage -python
"""N609X: normal test at the directly reconstructed H018 slope root."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import Matrix, PolynomialRing, vector
from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def lead(x): return x[x.valuation()]
def coeffs(field,pts,vals): return Matrix(field,[[field(s)**i for i in range(10)] for s in pts]).solve_right(vector(field,vals))

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
 target,cover,pi,phi=fixture(); field,ring,u,support,matrix,xvalues=origin_data(target,cover,pi,phi); det=matrix.det(); c=lead((-support[0]/support[1])/u); formal=target.formal_group(); rows=[]
 for level in (0,1,17):
  moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); samples=[]
  for n in range(1,32):
   slope=field(n)
   if slope==c: continue
   parameter=ring(slope)*u
   try: raw=sum(moving[i]*moving_basis((ring(formal.x(32)(parameter)),ring(formal.y(32)(parameter))),support)[i] for i in range(9))
   except ZeroDivisionError: continue
   normal=-det*parameter**8*raw
   if normal.valuation()>=0: samples.append((slope,normal[0]))
  cs=coeffs(field,[s for s,_ in samples[:10]],[(s-c)*v for s,v in samples[:10]]); R=PolynomialRing(field,'s'); s=R.gen(); num=sum(cs[i]*s**i for i in range(10)); roots=num.roots(field)
  for root,mult in roots:
   parameter=ring(root)*u; raw=sum(moving[i]*moving_basis((ring(formal.x(32)(parameter)),ring(formal.y(32)(parameter))),support)[i] for i in range(9)); normal=-det*parameter**8*raw
   rows.append({'level':level,'root':str(root),'root_multiplicity':int(mult),'normal_valuation':int(normal.valuation()),'normal_leading_coefficient':str(normal[normal.valuation()])})
 gates={'one_direct_root_per_level':len(rows)==3,'all_direct_roots_simple':all(row['root_multiplicity']==1 for row in rows),'all_direct_roots_have_uniform_order_four_normal_contact':all(row['normal_valuation']==4 for row in rows)}
 o={'schema':'ecdlp.h018.poincare-direct-slope-normal.n609x.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / DIRECT H018 WEIGHTED-BLOWUP RATIONAL STRICT-TRANSFORM SCREEN / MODEL-BOUND / TOY-EVIDENCE / RECIPROCAL_CHART_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'Each directly reconstructed exceptional slope root is simple, while the normal correction begins at uniform order four. The simple slope derivative screens the rational strict transform in this chart despite higher normal contact.','next_requirement':'Test reciprocal and geometric slope charts before global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}
 a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
 if not o['preflight_pass']: raise RuntimeError('N609X direct normal test failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()
