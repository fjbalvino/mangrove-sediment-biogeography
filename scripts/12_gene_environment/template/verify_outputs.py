"""Verify scientific invariants using the actual input release, never fixed n."""
from pathlib import Path
import argparse
import importlib.util
import json
import re
import numpy as np
import pandas as pd


def verify(root, viewers=False, cycle=None):
    src=root/'atlas/results'
    ids=(root/'inputs/functional_accessions.txt').read_text().splitlines()
    mat=pd.read_csv(root/'inputs/functional_KO_counts_original.tsv.gz',sep='\t',index_col=0).loc[:,ids]
    raw=pd.read_csv(src/'raw_assignments.tsv.gz',sep='\t',index_col=0)
    rel=pd.read_csv(src/'relative_representation.tsv.gz',sep='\t',index_col=0)
    assert len(ids)==len(set(ids)) and list(raw.columns)==ids
    assert raw.index.is_unique and rel.index.is_unique
    np.testing.assert_allclose(raw,mat.reindex(raw.index),rtol=0,atol=0,equal_nan=True)
    np.testing.assert_allclose(rel,mat.reindex(rel.index).div(mat.sum(),axis=1)*1e6,rtol=2e-13,atol=1e-9,equal_nan=True)
    missing=raw.index.difference(mat.index)
    assert raw.loc[missing].isna().all().all()
    weights=pd.read_csv(src/'unit_design_manifest.tsv',sep='\t',index_col=0)
    actual=pd.read_csv(src/'group_relative_representation.tsv.gz',sep='\t',index_col=0)
    assert weights.index.is_unique and set(weights.index)==set(ids)
    for group,rows in weights.groupby('analysis_group'):
        np.testing.assert_allclose(rows.descriptive_run_weight_within_group.sum(),1,atol=1e-12)
        expected=rel[rows.index].mul(rows.descriptive_run_weight_within_group,axis=1).sum(axis=1,min_count=1)
        np.testing.assert_allclose(actual.loc[group,expected.index],expected,rtol=2e-13,atol=1e-9,equal_nan=True)
    rules=pd.read_csv(root/'inputs/operational_rules.tsv',sep='\t')
    spec=importlib.util.spec_from_file_location('rules',root/'evaluate.py');ev=importlib.util.module_from_spec(spec);spec.loader.exec_module(ev)
    ev.self_tests()
    states=[]
    blocks=pd.read_csv(root/'outputs/rule_block_evidence.tsv.gz',sep='\t')
    assert not blocks.duplicated(['run_accession','rule_id','threshold']).any()
    assert len(blocks)==len(ids)*len(rules)*2
    assert set(blocks.run_accession)==set(ids) and set(blocks.rule_id)==set(rules.rule_id)
    assert set(blocks.threshold)=={1,2}
    assert (blocks.n_blocks>0).all()
    assert (blocks[['n_detected','n_not_met','n_unknown']]>=0).all().all()
    np.testing.assert_array_equal(blocks.n_detected+blocks.n_not_met+blocks.n_unknown,blocks.n_blocks)
    np.testing.assert_allclose(blocks.coverage_lower,blocks.n_detected/blocks.n_blocks)
    np.testing.assert_allclose(blocks.coverage_upper,(blocks.n_detected+blocks.n_unknown)/blocks.n_blocks)
    for threshold in [1,2]:
        block_states=pd.read_csv(root/f'outputs/rule_block_states_ge{threshold}.tsv.gz',sep='\t').set_index(['rule_id','block'])
        assert block_states.index.is_unique and list(block_states.columns)==ids
        state=pd.read_csv(src/f'rule_status_ge{threshold}.tsv',sep='\t',index_col=0).loc[ids]
        values={k:(mat.loc[k].to_numpy()>=threshold).astype('uint8') for k in raw.index if k in mat.index}
        for row in rules.itertuples():
            expected_blocks=ev.block_evidence(row.expression,values,len(ids))
            np.testing.assert_array_equal(block_states.loc[row.rule_id].to_numpy(),expected_blocks['states'])
            np.testing.assert_array_equal(state[row.rule_id],ev.evaluate(ev.parse(row.expression),values,len(ids)))
            b=blocks[(blocks.threshold==threshold)&(blocks.rule_id==row.rule_id)].set_index('run_accession').loc[ids]
            implied=np.where(b.n_not_met>0,0,np.where(b.n_unknown>0,2,1))
            np.testing.assert_array_equal(state[row.rule_id],implied)
        states.append(state.to_numpy())
    assert (states[0][states[1]==1]==1).all()
    diagnostics=pd.read_csv(src/'environment_association_diagnostics.tsv',sep='\t')
    assert not set(diagnostics)&{'p','q','p_value','q_value'}
    assert not diagnostics.duplicated(['KO','outcome','predictor']).any()
    cycles=json.loads((root/'atlas/curation/cycles.json').read_text())
    assert len(cycles)==8 and 'N' in cycles
    if viewers:
        key=pd.read_csv(src/'heatmap_sample_key.tsv',sep='\t')
        assert set(key.run_accession)==set(ids) and (np.diff(key.lon)<=0).all()
        assert list(key.column_left_to_right)==list(range(1,len(ids)+1))
        for c,info in cycles.items():
            if cycle and c!=cycle:continue
            text=(root/'atlas/figures'/info['slug']/'viewer.html').read_text()
            assert '<html lang="en">' in text and 'id="map"' in text
            assert 'id="sensitivity-image"' in text and 'id="heatmap-image"' in text
            assert 'id="terrain-retry"' in text and '"offline_relief"' in text
            assert not re.findall(r'__[A-Z_]+__',text)
    result={'passed':True,'profiles':len(ids),'rules':len(rules),'panels':len(cycles),'unavailable_KOs':len(missing),'viewers_checked':viewers,'viewer_cycle':cycle or 'all',
            'scope':'Numerical and structural verification; not validation of ecological inference or browser rendering'}
    (root/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--viewers',action='store_true')
    parser.add_argument('--cycle',choices=['N','C','CH4_C1','S','P','Fe','H2_electrons','O2_light'])
    args=parser.parse_args()
    verify(Path(__file__).resolve().parent,args.viewers,args.cycle)
