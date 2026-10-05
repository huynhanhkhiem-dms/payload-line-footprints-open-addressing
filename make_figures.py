#!/usr/bin/env python3
import csv, sys
from collections import defaultdict
import matplotlib.pyplot as plt
INPUT=sys.argv[1] if len(sys.argv)>1 else 'results_main.csv'
SCHEMES=[('up','RP'),('dh','DH'),('wtl','WTL-DH'),('eh','EH'),('ehb','EHB'),('fh','FH')]
with open(INPUT,newline='',encoding='utf-8') as f: rows=list(csv.DictReader(f))
by=defaultdict(list); allowed=dict(SCHEMES)
for r in rows:
    if r['scheme'] not in allowed: continue
    q=dict(r); q['delta']=float(r['delta'])
    for k in ('ins_probes','ins_lines','pos_lines','neg_lines'): q[k]=float(r[k])
    by[r['scheme']].append(q)
for k in by: by[k].sort(key=lambda r:r['delta'])
def draw(metric,ylabel,name):
    fig,ax=plt.subplots(figsize=(4.9,3.25))
    markers=['o','s','^','D','v','P']
    linestyles=['-','--','-.',':','-','--']
    for i,(s,lbl) in enumerate(SCHEMES):
        rs=by.get(s,[])
        if not rs: continue
        ax.plot([r['delta'] for r in rs],[r[metric] for r in rs],marker=markers[i],linestyle=linestyles[i],linewidth=1.2,markersize=3.4,label=lbl)
    ax.set_xscale('log'); ax.invert_xaxis(); ax.set_yscale('log')
    ax.set_xlabel(r'Empty fraction $\delta$ (occupancy increases $\rightarrow$)')
    ax.set_ylabel(ylabel); ax.grid(True,linewidth=0.35,alpha=0.35)
    ax.legend(ncol=3,fontsize=6.7,frameon=False); fig.tight_layout(); fig.savefig(name,bbox_inches='tight'); plt.close(fig)
draw('ins_probes','Mean slot probes / insertion','fig_build_probes.pdf')
draw('ins_lines','Mean payload lines / insertion','fig_build_lines.pdf')
draw('pos_lines','Mean payload lines / successful lookup','fig_positive_lines.pdf')
draw('neg_lines','Mean payload lines / unsuccessful lookup','fig_negative_lines.pdf')
