"""Real OP3 release acceptance: installed CLI, replay, portable inference, audits.

No model development or test-driven hyperparameter selection happens here.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from perturb_predict.io import digest, json_write, provenance


def call(*args):
    cp = subprocess.run([sys.executable, '-m', 'perturb_predict', *map(str,args)],
                        text=True, capture_output=True)
    if cp.returncode:
        raise RuntimeError(cp.stderr)
    return json.loads(cp.stdout)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    result={}
    for iteration in (1,2):
        root=a.out/f'run{iteration}'
        result[str(iteration)]=call('benchmark-op3','--data',a.data,'--out',root,
                                   '--accept-conditional','--accept-public-reuse')
    differences={}
    for name,n_cal,n_query in [('B_cells',11,49),('Myeloid_cells',10,47)]:
        r1,r2=a.out/'run1',a.out/'run2'
        with np.load(r1/name/'prediction/predictions.npz',allow_pickle=False) as z: p=z['prediction']
        with np.load(r2/name/'prediction/predictions.npz',allow_pickle=False) as z: q=z['prediction']
        np.testing.assert_array_equal(p,q)
        assert p.shape==(n_query,5288)
        fit=json.loads((r1/name/'model/fit.json').read_text())
        assert fit['n_calibration_compounds']==n_cal
        for run in (r1,r2):
            audit=json.loads((run/'prepared/audit.json').read_text())
            assert all(not (set(row['splits']) & {'private_test','public_test'}) for row in audit['reader_log'])
            score=json.loads((run/name/'evaluation/metrics.json').read_text())
            assert all('private_test' not in row['splits'] for row in score['reader_log'])
        # The subprocess receives only a portable model and an unlabeled query bundle.
        portable=a.out/'portable'/name
        call('predict','--model',r1/name/'model/model.npz',
             '--data',r1/'prepared'/name/'query.npz','--out',portable)
        with np.load(portable/'predictions.npz',allow_pickle=False) as z:
            np.testing.assert_array_equal(p,z['prediction'])
        differences[name]={'repeat_max_abs':float(np.max(abs(p-q))),
                           'portable_inference_max_abs':0.,'shape':list(p.shape)}
    # Reference belongs to an already-used public set: regression acceptance only.
    expected=0.1980665327273829
    assert abs(result['1']['model']['mse']-expected)<1e-9
    json_write(a.out/'acceptance.json',{'status':'passed','checks':differences,
          'scores':result['1'],'data_sha256':digest(a.data),
          'independent_biological_validation':False,
          'installed_package':__import__('perturb_predict').__file__,**provenance()})
    print(json.dumps(result['1'],indent=2))

if __name__=='__main__': main()
