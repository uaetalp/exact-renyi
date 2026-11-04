import json
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT),str(ROOT/'src')]


class AccountingTests(unittest.TestCase):
    def test_earlier_implementations_against_labelled_oracle(self):
        from labeled_path_oracle import frozen
        from workload_core import calc as workload_calc
        from confirmation_core import calc as confirmation_calc
        reference=frozen()
        for b in (2,3):
            for arm in ('U','E'):
                spec=dict(a=3,b=b,k=2,T=3,p=.125,sigma=2,arm=arm)
                expected=float(reference[str(b)]['d120']['logF'])
                for calc in (workload_calc,confirmation_calc):
                    self.assertAlmostEqual(calc(spec)['logF'],expected,places=12)

    def test_deterministic_gaussian_shift(self):
        from workload_core import calc
        for a in (2,4):
            for arm in ('U','E'):
                result=calc(dict(a=a,b=2,k=1,T=3,p=1.,sigma=2,arm=arm,identity_C=True))
                self.assertAlmostEqual(result['logF'],a*(a-1)/4,places=12)

    def test_interpolation_endpoint_attribution(self):
        from accounting_portfolio import interpolated,convert
        candidates=[c for c in interpolated({4:3.,8:14.},.001) if c['bracket']==[4,8]]
        low=next(c for c in candidates if c['evaluation_order']==4)
        high=next(c for c in candidates if c['evaluation_order']==8)
        middle=next(c for c in candidates if 4<c['evaluation_order']<8)
        self.assertEqual(low['active_orders'],[4])
        self.assertEqual(high['active_orders'],[8])
        self.assertEqual(low['epsilon'],convert(1.,4,.001))
        self.assertEqual(middle['active_orders'],[4,8])
        self.assertAlmostEqual(middle['evaluation_order'],math.exp(-5.25)/.001)

    def test_fixed_certificate_rejects_different_target_instance(self):
        from seed_core import certificate
        spec=json.loads((ROOT/'configs/main.json').read_text())['certificate_B']
        for extra in ({'a':2},{'T':4},{'C':'alternating_rows'}):
            with self.assertRaises(ValueError):
                certificate(dict(spec,**extra),'E')

    def test_recorded_results_and_paper_table_values(self):
        from summarize import tables
        data=tables()
        self.assertEqual(len(data['paired_cpu']),48)
        summary={r['case']:r for r in data['summary']}
        expected={'B':43.217,'L96':32.138,'L48':47.733,'BW4_LOW':21.484,'BW4_EQ':48.943,'P_HIGH':47.072,'C_VARY':46.835,'N2':-.078,'W0':.940,'W3':57.104,'W4':93.336}
        for case, value in expected.items():
            self.assertEqual(round(summary[case]['reduction_percent'],3),value)
        query=next(r for r in data['workload'] if r['case']=='W3' and float(r['delta'])==1e-5)
        self.assertEqual(round(query['B_epsilon'],6),4.179770)
        self.assertEqual(round(query['UE_epsilon'],6),.656104)

    def test_runner_order_matches_recorded_experiment_order(self):
        from reproduce import tasks_for
        main=json.loads((ROOT/'results/paper/main/calls.json').read_text())
        recorded=[c['id'] for c in main if c['group']!='A']
        generated=[]
        for suite in ('certified','expansion'):
            for t in tasks_for(suite,3):
                generated.append(f"B_R{t['round']}_{t['arm']}" if suite=='certified' else f"C_R{t['round']}_{t['case']}_{t['arm']}")
        self.assertEqual(generated,recorded)
        for suite in ('workload','confirmation'):
            filename='calls.json' if suite=='workload' else 'cost_summary.json'
            data=json.loads((ROOT/'results/paper'/suite/filename).read_text())
            calls=data if suite=='workload' else data['results']
            recorded=[tuple(c['name'].split('_')[i] for i in (0,2,3)) for c in calls if c['name'].startswith('R')]
            generated=[(f"R{t['round']}",t['case'],t['arm']) for t in tasks_for(suite,3)]
            self.assertEqual(generated,recorded)


if __name__=='__main__':
    unittest.main()
