#!/usr/bin/env sage -python
"""N609Y: reciprocal exceptional-direction screen from direct slope data."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import Matrix, PolynomialRing, vector
from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture
def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def lead(x): return x[x.valuation()]
def main():
 p=argparse.ArgumentParser(); p.add_argument('--out',type=Path,required=True); a=p.parse_args(); target,cover,pi,phi=fixture(); field,ring,u,support,matrix,xvalues=origin_data(target,cover,pi,phi); det=matrix.det(); c=lead((-support[0]/support[1])/u); formal=target.formal_group(); rows=[]
 for level in (0,1,17):
  moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); samples=[]
  for n in range(1,32):
   s=field(n)
   if s==c: continue
   q=ring(s)*u
   try: raw=sum(moving[i]*moving_basis((ring(formal.x(32)(q)),ring(formal.y(32)(q))),support)[i] for i in range(9))
   except ZeroDivisionError: continue
   value=-det*q**8*raw
   if value.valuation()>=0: samples.append((s,value[0]))
  coeff=Matrix(field,[[s**i for i in range(10)] for s,_ in samples[:10]]).solve_right(vector(field,[(s-c)*v for s,v in samples[:10]])); R=PolynomialRing(field,'s'); S=R.gen(); num=sum(coeff[i]*S**i for i in range(10)); rows.append({'level':level,'numerator_degree':int(num.degree()),'infinity_numerator_leading_coefficient':str(coeff[num.degree()]),'infinity_value_nonzero':coeff[num.degree()]!=0})
 gates={'three_levels_used':len(rows)==3,'all_direct_numerators_degree_one':all(r['numerator_degree']==1 for r in rows),'reciprocal_infinity_direction_is_nonzero':all(r['infinity_value_nonzero'] for r in rows)}
 o={'schema':'ecdlp.h018.poincare-direct-slope-infinity.n609y.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / RECIPROCAL EXCEPTIONAL DIRECTION SCREEN FROM DIRECT H018 SLOPE FORM / MODEL-BOUND / TOY-EVIDENCE / GEOMETRIC_ROOTS_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'The directly reconstructed slope form has degree-one numerator and a nonzero value at slope infinity, so the reciprocal exceptional direction is not a strict-transform root on the tested chart.','next_requirement':'Extend the direct slope reconstruction to geometric coefficients and prove chart transitions before global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}; a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii');
 if not o['preflight_pass']: raise RuntimeError('N609Y infinity screen failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()
