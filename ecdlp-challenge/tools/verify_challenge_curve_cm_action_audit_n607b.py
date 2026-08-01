#!/usr/bin/env sage -python
import argparse,json
from pathlib import Path
from math import isqrt
def main():
 p=argparse.ArgumentParser();p.add_argument("--primary",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();d=json.loads(a.primary.read_text());P=616883774851;T=984079;N=616882790773;m=P-(T*T-1)//4;c={"primary":d["preflight_pass"] is True,"norm":d["records"]["minimum_norm_for_b_nonzero"]==m,"rho":m>isqrt(N),"scalar_boundary":d["claim_status"].startswith("RESTRICTED THEOREM")};o={"verified":all(c.values()),"checks":c};a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps({"verified":o["verified"],"output":str(a.output)}))
if __name__=="__main__":main()
