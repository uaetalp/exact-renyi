#!/usr/bin/env python3
"""Run fixed paper configurations in fresh processes; no external dependencies."""
import argparse
from decimal import Decimal
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def config(suite):
    return read(ROOT / 'configs' / (suite + '.json'))


def compact(value):
    if isinstance(value, dict):
        return {k: compact(v) for k, v in value.items() if k not in {'trace', 'weights'}}
    if isinstance(value, list):
        return [compact(v) for v in value]
    return value


def worker(task, output, timeout):
    from runtime import configure
    configure(timeout)
    suite, case, arm = (task[x] for x in ('suite', 'case', 'arm'))
    trace = output / (task['id'] + '_traces')
    trace.mkdir()
    if suite == 'qa':
        from qa_seed import qa
        q = next(q for q in config('main')['A_QA'] if q['id'] == case)
        result = qa(q, trace)
    elif suite in ('certified', 'expansion'):
        from seed_core import calc, certificate
        cfg = config('main')
        if suite == 'certified':
            result = certificate(cfg['certificate_B'], arm)
        else:
            spec = next(s for s in cfg['C_cases'] if s['id'] == case)
            result = calc(spec, arm, trace=trace)
    elif suite == 'workload':
        from accounting_portfolio import workload
        cfg = config('workload')
        result = workload(cfg['configs'][case], arm, cfg['deltas'])
    else:
        from confirmation_core import calc
        result = calc(dict(config('confirmation')['configs'][case], arm=arm))
    write(output / (task['id'] + '.json'), result)


def near(x, y):
    return abs(x-y) <= 1e-10 + 1e-9*abs(y)


def compare_pair(suite, u, e):
    if suite == 'certified':
        for result in (u, e):
            assert result['certificate_status'] == 'PASS'
            assert Decimal(result['forward_delta']['hi']) <= Decimal('0.00001')
            assert Decimal(result['reverse_delta']['hi']) <= Decimal('0.00001')
        assert max(Decimal(u['value']['lo']), Decimal(e['value']['lo'])) <= min(Decimal(u['value']['hi']), Decimal(e['value']['hi']))
        return {'status': 'PASS', 'check': 'certificate_and_interval_overlap'}
    if suite == 'workload':
        for delta, q in u['queries'].items():
            assert near(q['epsilon'], e['queries'][delta]['epsilon'])
        for a, m in u['moments'].items():
            assert near(m['logF'], e['moments'][a]['logF'])
        return {'status': 'PASS', 'check': 'all_moments_and_privacy_queries'}
    assert near(u['value'] if suite == 'expansion' else u['logF'],
                e['value'] if suite == 'expansion' else e['logF'])
    assert len(u['rows']) == len(e['rows'])
    max_diff = 0.0
    for ur, er in zip(u['rows'], e['rows']):
        assert ur['row'] == er['row'] and ur['states'] == er['states']
        if suite == 'expansion':
            def states(row):
                return {tuple(tuple(z) for z in key): v for key, v in read(row['trace'])}
        else:
            def states(row):
                return {tuple(tuple(z) for z in item[:-1]): item[-1] for item in row['weights']}
        us, es = states(ur), states(er)
        assert us.keys() == es.keys()
        for key, value in us.items():
            assert near(value, es[key])
            max_diff = max(max_diff, abs(value-es[key]))
    return {'status': 'PASS', 'check': 'all_row_states', 'max_log_weight_difference': max_diff}


def tasks_for(suite, repetitions, case=None):
    if suite == 'qa':
        for q in config('main')['A_QA']:
            if case is None or q['id'] == case:
                yield dict(id=q['id'], suite=suite, case=q['id'], arm='', round=1)
        return
    cfg = config('main' if suite in ('certified', 'expansion') else suite)
    for rnd in range(1, repetitions+1):
        if suite == 'certified':
            cases = ['B']
            arms = cfg['certificate_B']['order'][rnd-1]
        elif suite == 'expansion':
            cases = cfg['C_task_order'][rnd-1]
            arms = cfg['C_arm_order_by_round'][rnd-1]
        else:
            cases = cfg['round' + str(rnd)]
            arms = ({1: ['B','U','E'], 2: ['E','U','B'], 3: ['U','B','E']} if suite == 'workload'
                    else {1: ['U','E'], 2: ['E','U'], 3: ['U','E']})[rnd]
        for position, name in enumerate(cases, 1):
            if case is not None and case != name:
                continue
            case_arms = (['U','E'] if position % 2 else ['E','U']) if suite == 'confirmation' and rnd == 3 else arms
            for arm in case_arms:
                yield dict(id=f'{suite}_R{rnd}_{name}_{arm}', suite=suite, case=name, arm=arm, round=rnd)


def execute(args):
    suites = ['certified','expansion','workload','confirmation'] if args.suite == 'all' else [args.suite]
    if args.suite == 'all' and args.case:
        raise ValueError('--case applies to a single suite')
    tasks = [task for suite in suites for task in tasks_for(suite, args.repetitions, args.case)]
    if not tasks:
        raise ValueError('unknown case for this suite')
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    state = dict(status='RUNNING', python=sys.version.split()[0], platform=sys.platform,
                 timing='complete_child_user_plus_system_CPU', calls=[], comparisons=[])
    pairs = {}
    try:
        for task in tasks:
            before = resource.getrusage(resource.RUSAGE_CHILDREN)
            start = time.monotonic()
            env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0')
            with (out / (task['id'] + '.stderr')).open('w') as err:
                run = subprocess.run([sys.executable, str(Path(__file__).resolve()), '_worker',
                                      json.dumps(task), str(out), str(args.timeout)],
                                     stdout=subprocess.DEVNULL, stderr=err, env=env,
                                     timeout=args.timeout+5, check=False)
            after = resource.getrusage(resource.RUSAGE_CHILDREN)
            receipt = dict(task, cpu=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                           wall=time.monotonic()-start, returncode=run.returncode)
            state['calls'].append(receipt)
            write(out/'run.json', state)
            if run.returncode:
                raise RuntimeError(f"{task['id']} failed; see its stderr file")
            result = read(out/(task['id']+'.json'))
            assert result.get('status','PASS') == 'PASS'
            if task['arm'] in ('U','E'):
                key = (task['suite'], task['case'], task['round'])
                pair = pairs.setdefault(key, {})
                pair[task['arm']] = (task, result)
                if len(pair) == 2:
                    comparison = compare_pair(task['suite'], pair['U'][1], pair['E'][1])
                    state['comparisons'].append(dict(suite=key[0], case=key[1], round=key[2], **comparison))
                    if not args.keep_traces:
                        for old_task, data in pair.values():
                            write(out/(old_task['id']+'.json'), compact(data))
                            shutil.rmtree(out/(old_task['id']+'_traces'))
                    del pairs[key]
            elif not args.keep_traces:
                write(out/(task['id']+'.json'), compact(result))
                shutil.rmtree(out/(task['id']+'_traces'))
            print(f"{task['id']}: PASS ({receipt['cpu']:.3f} CPU s)", flush=True)
            write(out/'run.json', state)
        state['status'] = 'PASS'
    except Exception as exc:
        state.update(status='FAILED', error=str(exc))
        raise
    finally:
        write(out/'run.json', state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('suite', choices=['qa','certified','expansion','workload','confirmation','all'])
    parser.add_argument('--case', help='one case, e.g. L96, W0 or N2')
    parser.add_argument('--repetitions', type=int, choices=[1,2,3], default=3)
    parser.add_argument('--output', required=True, help='new directory for this run')
    parser.add_argument('--timeout', type=int, default=180, help='seconds per worker (default: 180)')
    parser.add_argument('--keep-traces', action='store_true', help='retain full generated row states')
    args = parser.parse_args()
    if args.timeout < 1:
        parser.error('--timeout must be positive')
    execute(args)


if __name__ == '__main__':
    if sys.flags.optimize:
        raise SystemExit('Run without -O: numerical checks use assertions.')
    if len(sys.argv) > 1 and sys.argv[1] == '_worker':
        worker(json.loads(sys.argv[2]), Path(sys.argv[3]), int(sys.argv[4]))
    else:
        main()
