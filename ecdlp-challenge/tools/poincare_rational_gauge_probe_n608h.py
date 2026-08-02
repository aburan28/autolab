#!/usr/bin/env sage -python
"""N608H: bounded-pole rational fit for the Poincare cycle gauge."""
import argparse, datetime, hashlib, json, sys
from pathlib import Path
from sage.all import EllipticCurve, GF, Matrix, gcd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import product_kummer_h018_poincare_scalar_cocycle_n606z as transport


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def basis(point, d):
    x, y = point[0], point[1]
    return [x**a * (y if b else 1) for b in (0, 1) for a in range(d + 1) if 2*a + 3*b <= d]

def main():
    p=argparse.ArgumentParser(); p.add_argument("--out",type=Path,required=True); p.add_argument("--samples",type=int,default=24); a=p.parse_args()
    base=EllipticCurve(GF(103),[0,0,0,1,24]); field=GF(103**2,name="a"); curve=base.base_extend(field)
    g0 = next(q for q in base.points() if not q.is_zero())
    g = curve(field(g0[0]), field(g0[1]))
    orbit=[i*(84*g) for i in range(109)]; inverse=pow(109,-1,field.order()-1); data=[]
    for point in curve.points():
        if point.is_zero(): continue
        try: value=field.prod(transport.transport(orbit[i],g,point+(108-i)*g) for i in range(109))
        except ZeroDivisionError: continue
        if value:
            data.append((point,value**inverse))
        if len(data)==a.samples: break
    fits=[]
    for d in (2,4,6):
        width=len(basis(data[0][0],d)); rows=[]
        for point,value in data:
            values=basis(point,d); rows.append(values+[-value*x for x in values])
        kernel=Matrix(field,rows).right_kernel(); valid=False
        for v in kernel.basis():
            num,den=v[:width],v[width:]
            if any(num) and any(den) and all(sum(c*x for c,x in zip(den,basis(point,d))) != 0 for point,_ in data): valid=True
        fits.append({"pole_bound":d,"basis_dimension":width,"equation_count":len(rows),"nullity":int(kernel.dimension()),"has_regular_nonconstant_fit":valid})
    gates={"enough_regular_samples":len(data)==a.samples,"unique_roots":gcd(109,field.order()-1)==1,"no_low_pole_rational_fit":not any(row["has_regular_nonconstant_fit"] for row in fits)}
    out={"schema":"ecdlp.h018.poincare-rational-gauge.n608h.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha(Path(__file__)),"claim_status":"HYPOTHESIS / BOUNDED_POLE_POINCARE_GAUGE_PREFLIGHT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM","records":{"field_degree":2,"sample_count":len(data),"fits":fits},"gates":gates,"preflight_pass":all(gates.values()),"next_requirement":"A surviving rational gauge needs geometric divisor/normalization proof and pushforward-frame comparison; a negative only rejects these bounded pole spaces."}
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="ascii"); print(json.dumps({"checks":sum(gates.values()),"output":str(a.out)},sort_keys=True))

if __name__=="__main__": main()
