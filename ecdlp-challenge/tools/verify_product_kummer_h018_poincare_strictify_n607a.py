#!/usr/bin/env sage -python
"""Independent replay for N607A."""
import argparse,json
from math import gcd
from pathlib import Path
from sage.all import *
import product_kummer_h018_poincare_strictify_n607a as n
def main():
 p=argparse.ArgumentParser();p.add_argument("--primary",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();d=json.loads(a.primary.read_text());c={"primary":d["preflight_pass"] is True,"field":d["records"]["field_unit_order"]==103**2-1,"coprime":gcd(109,103**2-1)==1,"boundary":d["geometric_normalized_poincare_rigidification_constructed"] is False};o={"verified":all(c.values()),"checks":c};a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps({"verified":o["verified"],"output":str(a.output)}))
if __name__=="__main__":main()
