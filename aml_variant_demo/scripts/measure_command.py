#!/usr/bin/env python3
"""Run one command and report wall time and child-process maximum RSS on macOS."""
import argparse, json, resource, subprocess, time
ap=argparse.ArgumentParser(); ap.add_argument('--label',required=True); ap.add_argument('--output',required=True); ap.add_argument('command',nargs=argparse.REMAINDER); a=ap.parse_args()
if a.command and a.command[0]=='--': a.command=a.command[1:]
t=time.monotonic(); p=subprocess.run(a.command); wall=time.monotonic()-t
# Darwin reports ru_maxrss in bytes; Linux reports KiB.
rss=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
with open(a.output,'w') as f: json.dump({'label':a.label,'exit':p.returncode,'wall_seconds':round(wall,3),'peak_rss_bytes':rss},f,indent=2)
raise SystemExit(p.returncode)
