import json, sys, math
S = sys.argv[1]
for st in sys.argv[2:]:
    with open(f'{S}/grade/v5-{st}/results.json', encoding='utf-8') as result_file:
        results = json.load(result_file)
    for rr in results:
        if 'error' in rr: print(rr['id'], 'ERROR', rr['error']); continue
        o, re, d, rl = rr['ours'], rr['real'], rr['distance'], rr['rulers']
        f = lambda v: 'nan' if v is None or (isinstance(v, float) and math.isnan(v)) else round(v, 2)
        print(f"{rr['id']:24s} f0 {o['f0']:7.1f} dist {d['spec_db']:5.1f} (env {f(d['env_db'])} harm {f(d['harm_db'])}) "
              f"dynR {[round(v,1) for v in rl.get('other_dynamic',[])]} nbR {[round(v,1) for v in rl.get('neighbour',[])][:3]} "
              f"near {[(n[0], round(n[1],1)) for n in rr['nearest'][:3]]} HNR {round(o['hnr_db'],1)}/{round(re['hnr_db'],1)} "
              f"fluct {f(o.get('amp_fluct_db'))}/{f(re.get('amp_fluct_db'))} sub {f(o.get('subharmonic_db'))}/{f(re.get('subharmonic_db'))}")
