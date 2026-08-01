#!/usr/bin/env sage -python
import argparse,json
from pathlib import Path
def main():
 p=argparse.ArgumentParser();p.add_argument("--primary",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();d=json.loads(a.primary.read_text());c={"primary":d["preflight_pass"] is True,"transitions":d["records"]["transition_count"]==109,"invertible":d["records"]["monodromy_rank"]==9,"nonscalar":d["records"]["monodromy_is_scalar"] is False,"untwisted":d["theta11_twist_included"] is False};o={"verified":all(c.values()),"checks":c};a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps({"verified":o["verified"],"output":str(a.output)}))
if __name__=="__main__":main()
