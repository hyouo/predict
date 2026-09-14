"""Exploratory correlated-anchor fitting; see protocol 005."""
from __future__ import annotations
import argparse, json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from tools.assets import ROOT, sha256_file
from src.rna_data import OP3Store, TARGET_CONTEXTS, NAMES as OP3_NAMES
from src.rna_metrics import effect_metrics, row_correlation
from src.correlated_anchors import fit_correlated_bank
from tools.sciplex_experiment import Store as SciStore, CONTEXTS, NAMES as SCI_NAMES


def op3_factors(store, context, source, anchors, grouping):
    o = store.obs; donors = source.donors
    cc, _, _ = store._controls(context, donors, 'C', 'anchor_control')
    cd, _, _ = store._controls(context, donors, 'D', 'anchor_control')
    drugs = [source.drugs[i] for i in anchors.indices]
    rows = o[(o.cell_type == context) & o.donor.isin(donors) & (~o.is_control)
             & o['split'].eq('train') & o.perturbagen.isin(drugs)]
    if grouping not in ('donor', 'plate'):
        raise ValueError('Grouping must be donor or plate')
    units = sorted(rows[grouping].unique()); count = rows.groupby('perturbagen').size()
    f = np.zeros((len(units), len(drugs), source.means.shape[2]))
    for _, r in rows.iterrows():
        f[units.index(r[grouping]), drugs.index(r.perturbagen)] += (cd[r.plate]-cc[r.plate])/count[r.perturbagen]
    return f, units


def sciplex_factors(store, genes):
    controls, _, _ = store.controls(store.target, 'A', 'anchor_control', genes)
    o = store.o
    rows = o[o.qualifying & o.cell_type.eq(store.target) & o.panel.eq('A') & o.cid.isin(store.anchors)]
    units = sorted(rows.plate.unique()); count = rows.groupby('cid').size()
    f = np.zeros((len(units), len(store.anchors), len(genes)))
    for _, r in rows.iterrows():
        f[units.index(r.plate), store.anchors.index(r.cid)] += (controls[r.plate][1]-controls[r.plate][0])/count[r.cid]
    return f, units


def run(args):
    out = args.output; out.mkdir(parents=True, exist_ok=False)
    base = args.frozen; records = {}; contexts = {}; result = {}; sample_access = {}
    if args.dataset == 'op3':
        oldlock = json.loads((base/'prediction_lock.json').read_text())
        store = OP3Store(ROOT/'data/rna/op3_standardized_processed.h5ad',
                         json.loads((ROOT/'reports/data/op3_donors.json').read_text()))
        source = store.source(tuple(oldlock['donors_for_fitting']))
        for c in TARGET_CONTEXTS:
            name = OP3_NAMES[c]; a = store.anchors(c, source)
            original = json.loads((base/f'{name}_fit.json').read_text())['few_shot']
            result[name] = {}; records[name] = {}
            for mode in ('full', 'diagonal'):
                grouping = 'donor'
                factors, units = op3_factors(store, c, source, a, grouping)
                pred, rec = fit_correlated_bank(source.means, a.indices, a.reference_c, a.reference_d,
                                               factors, noise_mode=mode)
                result[name].update({mode+'__'+k: v for k, v in pred.items()})
                records[name][mode] = {'reference_units': units, **rec}
            contexts[name] = {'labels': source.drugs, 'anchor_count': len(a.indices)}
        sample_access['fit'] = store.access
        truth_dir = base/args.evaluation
        labels_key = 'compound'
    else:
        for c in CONTEXTS:
            name = SCI_NAMES[c]; sub = base/name; oldlock = json.loads((sub/'prediction_lock.json').read_text())
            store = SciStore(ROOT/'data/rna/srivatsan20_sciplex3_processed.h5ad', c)
            genes = np.asarray(oldlock['gene_indices']); inputs = np.load(sub/'inputs.npz', allow_pickle=False)
            original = json.loads((sub/'fit.json').read_text())['few_shot']
            factors, units = sciplex_factors(store, genes)
            result[name] = {}; records[name] = {}
            for mode in ('full', 'diagonal'):
                pred, rec = fit_correlated_bank(inputs['source'], inputs['anchor_indices'], inputs['anchor_c'],
                                                inputs['anchor_d'], factors, noise_mode=mode)
                result[name].update({mode+'__'+k: v for k, v in pred.items()})
                records[name][mode] = {'reference_units': units, **rec}
            sample_access[name] = store.access
            contexts[name] = {'labels': oldlock['cids'], 'anchor_count': 16}
        truth_dir = base/'evaluation'; labels_key = 'cid'
    hashes = {}
    for name, pred in result.items():
        file = out/f'{name}_predictions.npz'; np.savez_compressed(file, **pred); hashes[file.name] = sha256_file(file)
    lock = {'scientific_status': 'exploratory_followup_not_fresh_confirmation', 'dataset': args.dataset,
            'frozen_input_directory': str(base), 'prediction_hashes': hashes, 'contexts': contexts,
            'code_sha256': {p: sha256_file(ROOT/p) for p in ['src/correlated_anchors.py','src/reference_risk.py','tools/correlated_anchor_experiment.py','src/rna_models.py']},
            'predictions_sealed_before_this_evaluation': True, 'sample_access': sample_access,
            'sealed_utc': datetime.now(timezone.utc).isoformat()}
    (out/'prediction_lock.json').write_text(json.dumps(lock, indent=2)+'\n')
    (out/'fit.json').write_text(json.dumps(records, indent=2)+'\n')
    summaries = []; rows = []
    for name, predictions in result.items():
        truth = np.load(truth_dir/f'{name}_truth.npz', allow_pickle=False)
        t = truth['mean'] if args.dataset == 'op3' else truth['truth']
        ix = truth['source_indices'] if args.dataset == 'op3' else truth['query_indices']
        labels = contexts[name]['labels']
        for method, p in predictions.items():
            p = p[ix]; summaries.append({'context': name, 'method': method, **effect_metrics(p,t)})
            cor = row_correlation(p,t)
            for j, q in enumerate(ix):
                rows.append({'context': name, 'method': method, labels_key: labels[q],
                             'mse': float(np.mean((p[j]-t[j])**2)), 'pearson': float(cor[j])})
    d = pd.DataFrame(summaries); d.to_csv(out/'by_context.csv', index=False)
    pd.DataFrame(rows).to_csv(out/'by_compound.csv', index=False)
    d.groupby('method').mean(numeric_only=True).to_csv(out/'macro_summary.csv')
    print(d[['context','method','mse','centered_mse','bias_mse']].to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', choices=['op3','sciplex'], required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--evaluation', default='evaluation_private_test')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
