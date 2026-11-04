#!/usr/bin/env python3
"""Rebuild paper tables from the included recorded results, without rerunning experiments."""
import argparse
import csv
from decimal import Decimal
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(Path(path).read_text())


def near(x, y):
    return abs(x-y) <= 1e-10 + 1e-9*abs(y)


def tables(root=ROOT / 'results/paper'):
    rows = []
    main = root / 'main'
    calls = read(main/'calls.json')
    for c in calls:
        assert c['status'] == 'PASS' and c['returncode'] == 0
    assert len(calls) == 45
    costs = {c['id']: c['cpu'] for c in calls}
    for case in ['B','L96','L48','BW4_LOW','BW4_EQ','P_HIGH','C_VARY']:
        for rnd in (1,2,3):
            ids = {arm: f'B_R{rnd}_{arm}' if case == 'B' else f'C_R{rnd}_{case}_{arm}' for arm in ('U','E')}
            results = {arm: read(main/(ident+'.json')) for arm, ident in ids.items()}
            if case == 'B':
                for v in results.values():
                    assert v['certificate_status'] == 'PASS'
                    assert max(Decimal(v[x]['hi']) for x in ('forward_delta','reverse_delta')) <= Decimal('0.00001')
                    assert Decimal(v['rdp_width']) <= Decimal('1e-8')
                assert max(Decimal(v['value']['lo']) for v in results.values()) <= min(Decimal(v['value']['hi']) for v in results.values())
            else:
                assert near(results['U']['value'], results['E']['value'])
                comparison = read(main/f'{case}_R{rnd}_comparison.json')
                assert len(comparison) == len(results['U']['rows']) == len(results['E']['rows'])
                assert all(x['max_log_weight_difference'] <= 1e-9 for x in comparison)
            u, e = (costs[ids[arm]] for arm in ('U','E'))
            rows.append(dict(suite='certified' if case=='B' else 'expansion',case=case,pair=rnd,U_cpu_s=u,E_cpu_s=e,reduction_percent=100*(1-e/u)))

    work = root/'workload'
    wcalls = read(work/'calls.json')
    wcosts = {c['name']: c['cpu'] for c in wcalls}
    queries = []
    for case in ['W0','W3','W4']:
        for rnd in (1,2,3):
            files = {arm: next(work.glob(f'R{rnd}_*_{case}_{arm}.json')) for arm in ('B','U','E')}
            values = {arm:read(p) for arm,p in files.items()}
            assert all(v['status'] == 'PASS' for v in values.values())
            b, u, e = (wcosts[files[arm].stem] for arm in ('B','U','E'))
            rows.append(dict(suite='workload',case=case,pair=rnd,U_cpu_s=u,E_cpu_s=e,reduction_percent=100*(1-e/u)))
            for delta, q in values['U']['queries'].items():
                assert near(q['epsilon'], values['E']['queries'][delta]['epsilon'])
                queries.append(dict(case=case,pair=rnd,delta=delta,B_epsilon=values['B']['queries'][delta]['epsilon'],UE_epsilon=q['epsilon'],B_cpu_s=b,U_cpu_s=u,E_cpu_s=e))

    confirm = read(root/'confirmation/cost_summary.json')['results']
    assert len(confirm) == 36 and all(c['returncode'] == 0 for c in confirm)
    for case in ['N2','Z2','P2','N3','P3','A4']:
        for rnd in (1,2,3):
            pair = {c['arm']:c for c in confirm if c['config']==case and c['name'].startswith(f'R{rnd}_')}
            assert near(pair['U']['logF'], pair['E']['logF'])
            u,e = (pair[arm]['cpu'] for arm in ('U','E'))
            rows.append(dict(suite='confirmation',case=case,pair=rnd,U_cpu_s=u,E_cpu_s=e,reduction_percent=100*(1-e/u)))

    medians = []
    for suite, case in dict.fromkeys((r['suite'],r['case']) for r in rows):
        subset = [r for r in rows if (r['suite'],r['case']) == (suite,case)]
        medians.append(dict(suite=suite,case=case,U_cpu_s=statistics.median(r['U_cpu_s'] for r in subset),E_cpu_s=statistics.median(r['E_cpu_s'] for r in subset),reduction_percent=statistics.median(r['reduction_percent'] for r in subset)))
    workload = []
    for case, delta in dict.fromkeys((q['case'],q['delta']) for q in queries):
        subset = [q for q in queries if (q['case'],q['delta'])==(case,delta)]
        workload.append(dict(case=case,delta=delta,**{k:statistics.median(q[k] for q in subset) for k in ['B_epsilon','UE_epsilon','B_cpu_s','U_cpu_s','E_cpu_s']}))
    return {'paired_cpu':rows,'summary':medians,'privacy_queries':queries,'workload':workload}


def save(destination):
    data = tables()
    destination.mkdir(parents=True,exist_ok=True)
    for name, rows in data.items():
        with (destination/(name+'.csv')).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/derived')
    args=parser.parse_args()
    if __import__('sys').flags.optimize:
        raise SystemExit('Run without -O: result checks use assertions.')
    data=save(args.output)
    for row in data['summary']:
        print(f"{row['suite']:12s} {row['case']:8s} U={row['U_cpu_s']:.6f} E={row['E_cpu_s']:.6f} reduction={row['reduction_percent']:.3f}%")
